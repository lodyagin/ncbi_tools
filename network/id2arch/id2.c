#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <id2map.h>
#include <id2gen.h>

static Boolean loaded = FALSE;

#include <id2.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
id2genAsnLoad(void)
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
*    Generated object loaders for Module NCBI-ID2Access
*    Generated using ASNCODE Revision: 6.0 at Dec 15, 2003  5:08 PM
*
**************************************************/


/**************************************************
*
*    ID2RequestPacketFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestPacketPtr LIBCALL
ID2RequestPacketFree(ID2RequestPacketPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) ID2RequestFree);
   return NULL;
}


/**************************************************
*
*    ID2RequestPacketAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestPacketPtr LIBCALL
ID2RequestPacketAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestPacketPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2RequestPacket ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_PACKET);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_PACKET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2RequestAsnRead, (AsnOptFreeFunc) ID2RequestFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestPacketFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2RequestPacketAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestPacketAsnWrite(ID2RequestPacketPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_PACKET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) ID2RequestAsnWrite, aip, atp, ID2_REQUEST_PACKET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2RequestNew()
*
**************************************************/
NLM_EXTERN 
ID2RequestPtr LIBCALL
ID2RequestNew(void)
{
   ID2RequestPtr ptr = MemNew((size_t) sizeof(ID2Request));

   return ptr;

}


/**************************************************
*
*    ID2RequestFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestPtr LIBCALL
ID2RequestFree(ID2RequestPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2ParamsFree(ptr -> params);
   Request_requestFree(ptr -> Request_request);
   return MemFree(ptr);
}


/**************************************************
*
*    Request_requestFree()
*
**************************************************/
static 
Request_requestPtr LIBCALL
Request_requestFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case Request_request_get_packages:
      ID2RequestGetPackagesFree(anp -> data.ptrvalue);
      break;
   case Request_request_string_to_gi:
      ID2RequestStringToGiFree(anp -> data.ptrvalue);
      break;
   case Request_request_seq_id_to_gi:
      ID2RequestSeqIdToGiFree(anp -> data.ptrvalue);
      break;
   case Request_request_gi_to_tse_id:
      ID2RequestGiToTSEIdFree(anp -> data.ptrvalue);
      break;
   case Request_request_get_tse:
      ID2RequestGetTSEFree(anp -> data.ptrvalue);
      break;
   case Request_request_reget_tse:
      ID2RequestReGetTSEFree(anp -> data.ptrvalue);
      break;
   case Request_request_get_chunks:
      ID2SRequestGetChunksFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ID2RequestAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestPtr LIBCALL
ID2RequestAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2Request ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2RequestNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REQUEST_serial_number) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> serial_number = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REQUEST_params) {
      ptr -> params = ID2ParamsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REQUEST_request) {
      ptr -> Request_request = Request_requestAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestFree(ptr);
   goto ret;
}



/**************************************************
*
*    Request_requestAsnRead()
*
**************************************************/
static 
Request_requestPtr LIBCALL
Request_requestAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Request_request ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_request);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_request);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == ID2_REQUEST_request_init) {
      choice = Request_request_init;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == REQUEST_request_get_packages) {
      choice = Request_request_get_packages;
      func = (AsnReadFunc) ID2RequestGetPackagesAsnRead;
   }
   else if (atp == REQUEST_request_string_to_gi) {
      choice = Request_request_string_to_gi;
      func = (AsnReadFunc) ID2RequestStringToGiAsnRead;
   }
   else if (atp == REQUEST_request_seq_id_to_gi) {
      choice = Request_request_seq_id_to_gi;
      func = (AsnReadFunc) ID2RequestSeqIdToGiAsnRead;
   }
   else if (atp == REQUEST_request_gi_to_tse_id) {
      choice = Request_request_gi_to_tse_id;
      func = (AsnReadFunc) ID2RequestGiToTSEIdAsnRead;
   }
   else if (atp == ID2_REQUEST_request_get_tse) {
      choice = Request_request_get_tse;
      func = (AsnReadFunc) ID2RequestGetTSEAsnRead;
   }
   else if (atp == ID2_REQUEST_request_reget_tse) {
      choice = Request_request_reget_tse;
      func = (AsnReadFunc) ID2RequestReGetTSEAsnRead;
   }
   else if (atp == ID2_REQUEST_request_get_chunks) {
      choice = Request_request_get_chunks;
      func = (AsnReadFunc) ID2SRequestGetChunksAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ID2RequestAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestAsnWrite(ID2RequestPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> serial_number;
   retval = AsnWrite(aip, ID2_REQUEST_serial_number,  &av);
   if (ptr -> params != NULL) {
      if ( ! ID2ParamsAsnWrite(ptr -> params, aip, ID2_REQUEST_params)) {
         goto erret;
      }
   }
   if (ptr -> Request_request != NULL) {
      if ( ! Request_requestAsnWrite(ptr -> Request_request, aip, ID2_REQUEST_request)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Request_requestAsnWrite()
*
**************************************************/
static Boolean LIBCALL 
Request_requestAsnWrite(Request_requestPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ID2_REQUEST_request);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case Request_request_init:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, ID2_REQUEST_request_init, &av);
      break;
   case Request_request_get_packages:
      writetype = REQUEST_request_get_packages;
      func = (AsnWriteFunc) ID2RequestGetPackagesAsnWrite;
      break;
   case Request_request_string_to_gi:
      writetype = REQUEST_request_string_to_gi;
      func = (AsnWriteFunc) ID2RequestStringToGiAsnWrite;
      break;
   case Request_request_seq_id_to_gi:
      writetype = REQUEST_request_seq_id_to_gi;
      func = (AsnWriteFunc) ID2RequestSeqIdToGiAsnWrite;
      break;
   case Request_request_gi_to_tse_id:
      writetype = REQUEST_request_gi_to_tse_id;
      func = (AsnWriteFunc) ID2RequestGiToTSEIdAsnWrite;
      break;
   case Request_request_get_tse:
      writetype = ID2_REQUEST_request_get_tse;
      func = (AsnWriteFunc) ID2RequestGetTSEAsnWrite;
      break;
   case Request_request_reget_tse:
      writetype = ID2_REQUEST_request_reget_tse;
      func = (AsnWriteFunc) ID2RequestReGetTSEAsnWrite;
      break;
   case Request_request_get_chunks:
      writetype = ID2_REQUEST_request_get_chunks;
      func = (AsnWriteFunc) ID2SRequestGetChunksAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ID2ParamsFree()
*
**************************************************/
NLM_EXTERN 
ID2ParamsPtr LIBCALL
ID2ParamsFree(ID2ParamsPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) ID2ParamFree);
   return NULL;
}


/**************************************************
*
*    ID2ParamsAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ParamsPtr LIBCALL
ID2ParamsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ParamsPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2Params ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_PARAMS);
   } else {
      atp = AsnLinkType(orig, ID2_PARAMS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2ParamAsnRead, (AsnOptFreeFunc) ID2ParamFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ParamsFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ParamsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ParamsAsnWrite(ID2ParamsPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_PARAMS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) ID2ParamAsnWrite, aip, atp, ID2_PARAMS_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2RequestGetPackagesNew()
*
**************************************************/
NLM_EXTERN 
ID2RequestGetPackagesPtr LIBCALL
ID2RequestGetPackagesNew(void)
{
   ID2RequestGetPackagesPtr ptr = MemNew((size_t) sizeof(ID2RequestGetPackages));

   return ptr;

}


/**************************************************
*
*    ID2RequestGetPackagesFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestGetPackagesPtr LIBCALL
ID2RequestGetPackagesFree(ID2RequestGetPackagesPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> names ,ASNCODE_PTRVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2RequestGetPackagesAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestGetPackagesPtr LIBCALL
ID2RequestGetPackagesAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestGetPackagesPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2RequestGetPackages ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_GET_PACKAGES);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_GET_PACKAGES);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2RequestGetPackagesNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REQUEST_GET_PACKAGES_names) {
      ptr -> names = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> names == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GET_PACKAGES_no_contents) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> no_contents = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestGetPackagesFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2RequestGetPackagesAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestGetPackagesAsnWrite(ID2RequestGetPackagesPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_GET_PACKAGES);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> names ,ASNCODE_PTRVAL_SLOT, aip, ID2_REQUEST_GET_PACKAGES_names, REQUEST_GET_PACKAGES_names_E);
   av.boolvalue = ptr -> no_contents;
   retval = AsnWrite(aip, GET_PACKAGES_no_contents,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2RequestStringToGiNew()
*
**************************************************/
NLM_EXTERN 
ID2RequestStringToGiPtr LIBCALL
ID2RequestStringToGiNew(void)
{
   ID2RequestStringToGiPtr ptr = MemNew((size_t) sizeof(ID2RequestStringToGi));

   return ptr;

}


/**************************************************
*
*    ID2RequestStringToGiFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestStringToGiPtr LIBCALL
ID2RequestStringToGiFree(ID2RequestStringToGiPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> id);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2RequestStringToGiAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestStringToGiPtr LIBCALL
ID2RequestStringToGiAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestStringToGiPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2RequestStringToGi ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_STRING_TO_GI);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_STRING_TO_GI);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2RequestStringToGiNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REQUEST_STRING_TO_GI_id) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> id = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestStringToGiFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2RequestStringToGiAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestStringToGiAsnWrite(ID2RequestStringToGiPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_STRING_TO_GI);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> id != NULL) {
      av.ptrvalue = ptr -> id;
      retval = AsnWrite(aip, ID2_REQUEST_STRING_TO_GI_id,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2RequestSeqIdToGiNew()
*
**************************************************/
NLM_EXTERN 
ID2RequestSeqIdToGiPtr LIBCALL
ID2RequestSeqIdToGiNew(void)
{
   ID2RequestSeqIdToGiPtr ptr = MemNew((size_t) sizeof(ID2RequestSeqIdToGi));

   return ptr;

}


/**************************************************
*
*    ID2RequestSeqIdToGiFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestSeqIdToGiPtr LIBCALL
ID2RequestSeqIdToGiFree(ID2RequestSeqIdToGiPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SeqIdFree(ptr -> seq_id);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2RequestSeqIdToGiAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestSeqIdToGiPtr LIBCALL
ID2RequestSeqIdToGiAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestSeqIdToGiPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2RequestSeqIdToGi ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_SEQ_ID_TO_GI);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_SEQ_ID_TO_GI);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2RequestSeqIdToGiNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REQUEST_SEQ_ID_TO_GI_seq_id) {
      ptr -> seq_id = SeqIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestSeqIdToGiFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2RequestSeqIdToGiAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestSeqIdToGiAsnWrite(ID2RequestSeqIdToGiPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_SEQ_ID_TO_GI);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> seq_id != NULL) {
      if ( ! SeqIdAsnWrite(ptr -> seq_id, aip, ID2_REQUEST_SEQ_ID_TO_GI_seq_id)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2RequestGiToTSEIdNew()
*
**************************************************/
NLM_EXTERN 
ID2RequestGiToTSEIdPtr LIBCALL
ID2RequestGiToTSEIdNew(void)
{
   ID2RequestGiToTSEIdPtr ptr = MemNew((size_t) sizeof(ID2RequestGiToTSEId));

   return ptr;

}


/**************************************************
*
*    ID2RequestGiToTSEIdFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestGiToTSEIdPtr LIBCALL
ID2RequestGiToTSEIdFree(ID2RequestGiToTSEIdPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   Gi_giFree(ptr -> Gi_gi);
   AsnGenericBaseSeqOfFree(ptr -> sources ,ASNCODE_PTRVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    Gi_giFree()
*
**************************************************/
static 
Gi_giPtr LIBCALL
Gi_giFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case Gi_gi_string:
      ID2RequestStringToGiFree(anp -> data.ptrvalue);
      break;
   case Gi_gi_seq_id:
      ID2RequestSeqIdToGiFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ID2RequestGiToTSEIdAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestGiToTSEIdPtr LIBCALL
ID2RequestGiToTSEIdAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestGiToTSEIdPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2RequestGiToTSEId ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_GI_TO_TSE_ID);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_GI_TO_TSE_ID);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2RequestGiToTSEIdNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REQUEST_GI_TO_TSE_ID_gi) {
      ptr -> Gi_gi = Gi_giAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REQUEST_GI_TO_TSE_ID_sources) {
      ptr -> sources = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> sources == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == REQUEST_GI_TO_TSE_ID_external) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> external = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GI_TO_TSE_ID_current_gis) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> current_gis = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestGiToTSEIdFree(ptr);
   goto ret;
}



/**************************************************
*
*    Gi_giAsnRead()
*
**************************************************/
static 
Gi_giPtr LIBCALL
Gi_giAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Gi_gi ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_GI_TO_TSE_ID_gi);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_GI_TO_TSE_ID_gi);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == ID2_REQUEST_GI_TO_TSE_ID_gi_gi) {
      choice = Gi_gi_gi;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == REQUEST_GI_TO_TSE_ID_gi_string) {
      choice = Gi_gi_string;
      func = (AsnReadFunc) ID2RequestStringToGiAsnRead;
   }
   else if (atp == REQUEST_GI_TO_TSE_ID_gi_seq_id) {
      choice = Gi_gi_seq_id;
      func = (AsnReadFunc) ID2RequestSeqIdToGiAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ID2RequestGiToTSEIdAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestGiToTSEIdAsnWrite(ID2RequestGiToTSEIdPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_GI_TO_TSE_ID);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> Gi_gi != NULL) {
      if ( ! Gi_giAsnWrite(ptr -> Gi_gi, aip, ID2_REQUEST_GI_TO_TSE_ID_gi)) {
         goto erret;
      }
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> sources ,ASNCODE_PTRVAL_SLOT, aip, REQUEST_GI_TO_TSE_ID_sources, REQUEST_GI_TO_TSE_ID_sources_E);
   av.boolvalue = ptr -> external;
   retval = AsnWrite(aip, REQUEST_GI_TO_TSE_ID_external,  &av);
   av.boolvalue = ptr -> current_gis;
   retval = AsnWrite(aip, GI_TO_TSE_ID_current_gis,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Gi_giAsnWrite()
*
**************************************************/
static Boolean LIBCALL 
Gi_giAsnWrite(Gi_giPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ID2_REQUEST_GI_TO_TSE_ID_gi);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case Gi_gi_gi:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, ID2_REQUEST_GI_TO_TSE_ID_gi_gi, &av);
      break;
   case Gi_gi_string:
      writetype = REQUEST_GI_TO_TSE_ID_gi_string;
      func = (AsnWriteFunc) ID2RequestStringToGiAsnWrite;
      break;
   case Gi_gi_seq_id:
      writetype = REQUEST_GI_TO_TSE_ID_gi_seq_id;
      func = (AsnWriteFunc) ID2RequestSeqIdToGiAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ID2RequestGetTSENew()
*
**************************************************/
NLM_EXTERN 
ID2RequestGetTSEPtr LIBCALL
ID2RequestGetTSENew(void)
{
   ID2RequestGetTSEPtr ptr = MemNew((size_t) sizeof(ID2RequestGetTSE));

   return ptr;

}


/**************************************************
*
*    TseId_giNew()
*
**************************************************/
static 
TseId_giPtr LIBCALL
TseId_giNew(void)
{
   TseId_giPtr ptr = MemNew((size_t) sizeof(TseId_gi));

   return ptr;

}


/**************************************************
*
*    ID2RequestGetTSEFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestGetTSEPtr LIBCALL
ID2RequestGetTSEFree(ID2RequestGetTSEPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   TseId_tse_idFree(ptr -> TseId_tse_id);
   ID2GetTSEDetailsFree(ptr -> details);
   return MemFree(ptr);
}


/**************************************************
*
*    TseId_tse_idFree()
*
**************************************************/
static 
TseId_tse_idPtr LIBCALL
TseId_tse_idFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case TseId_tse_id_tse_id:
      ID2TSEIdFree(anp -> data.ptrvalue);
      break;
   case TseId_tse_id_TseId_Gi:
      TseId_giFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    TseId_giFree()
*
**************************************************/
static 
TseId_giPtr LIBCALL
TseId_giFree(TseId_giPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2RequestGiToTSEIdFree(ptr -> request);
   AsnGenericUserSeqOfFree(ptr -> exclude_tses, (AsnOptFreeFunc) ID2TSEIdFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2RequestGetTSEAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestGetTSEPtr LIBCALL
ID2RequestGetTSEAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestGetTSEPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2RequestGetTSE ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_GET_TSE);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_GET_TSE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2RequestGetTSENew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REQUEST_GET_TSE_tse_id) {
      ptr -> TseId_tse_id = TseId_tse_idAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REQUEST_GET_TSE_details) {
      ptr -> details = ID2GetTSEDetailsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestGetTSEFree(ptr);
   goto ret;
}



/**************************************************
*
*    TseId_tse_idAsnRead()
*
**************************************************/
static 
TseId_tse_idPtr LIBCALL
TseId_tse_idAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TseId_tse_id ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_GET_TSE_tse_id);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_GET_TSE_tse_id);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == REQUEST_GET_TSE_tse_id_tse_id) {
      choice = TseId_tse_id_tse_id;
      func = (AsnReadFunc) ID2TSEIdAsnRead;
   }
   else if (atp == ID2_REQUEST_GET_TSE_tse_id_gi) {
      choice = TseId_tse_id_TseId_Gi;
      func = (AsnReadFunc) TseId_giAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    TseId_giAsnRead()
*
**************************************************/
static 
TseId_giPtr LIBCALL
TseId_giAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TseId_giPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TseId_gi ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_GET_TSE_tse_id_gi);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_GET_TSE_tse_id_gi);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = TseId_giNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GET_TSE_tse_id_gi_request) {
      ptr -> request = ID2RequestGiToTSEIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GET_TSE_tse_id_gi_exclude_tses) {
      ptr -> exclude_tses = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2TSEIdAsnRead, (AsnOptFreeFunc) ID2TSEIdFree);
      if (isError && ptr -> exclude_tses == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = TseId_giFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2RequestGetTSEAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestGetTSEAsnWrite(ID2RequestGetTSEPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_GET_TSE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> TseId_tse_id != NULL) {
      if ( ! TseId_tse_idAsnWrite(ptr -> TseId_tse_id, aip, ID2_REQUEST_GET_TSE_tse_id)) {
         goto erret;
      }
   }
   if (ptr -> details != NULL) {
      if ( ! ID2GetTSEDetailsAsnWrite(ptr -> details, aip, ID2_REQUEST_GET_TSE_details)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    TseId_tse_idAsnWrite()
*
**************************************************/
static Boolean LIBCALL 
TseId_tse_idAsnWrite(TseId_tse_idPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ID2_REQUEST_GET_TSE_tse_id);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case TseId_tse_id_tse_id:
      writetype = REQUEST_GET_TSE_tse_id_tse_id;
      func = (AsnWriteFunc) ID2TSEIdAsnWrite;
      break;
   case TseId_tse_id_TseId_Gi:
      writetype = ID2_REQUEST_GET_TSE_tse_id_gi;
      func = (AsnWriteFunc) TseId_giAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    TseId_giAsnWrite()
*
**************************************************/
static Boolean LIBCALL 
TseId_giAsnWrite(TseId_giPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_GET_TSE_tse_id_gi);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> request != NULL) {
      if ( ! ID2RequestGiToTSEIdAsnWrite(ptr -> request, aip, GET_TSE_tse_id_gi_request)) {
         goto erret;
      }
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> exclude_tses, (AsnWriteFunc) ID2TSEIdAsnWrite, aip, GET_TSE_tse_id_gi_exclude_tses, TSE_tse_id_gi_exclude_tses_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2RequestReGetTSENew()
*
**************************************************/
NLM_EXTERN 
ID2RequestReGetTSEPtr LIBCALL
ID2RequestReGetTSENew(void)
{
   ID2RequestReGetTSEPtr ptr = MemNew((size_t) sizeof(ID2RequestReGetTSE));

   return ptr;

}


/**************************************************
*
*    ID2RequestReGetTSEFree()
*
**************************************************/
NLM_EXTERN 
ID2RequestReGetTSEPtr LIBCALL
ID2RequestReGetTSEFree(ID2RequestReGetTSEPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2TSEIdFree(ptr -> tse_id);
   ID2GetTSEDetailsFree(ptr -> details);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2RequestReGetTSEAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2RequestReGetTSEPtr LIBCALL
ID2RequestReGetTSEAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2RequestReGetTSEPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2RequestReGetTSE ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REQUEST_REGET_TSE);
   } else {
      atp = AsnLinkType(orig, ID2_REQUEST_REGET_TSE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2RequestReGetTSENew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REQUEST_REGET_TSE_tse_id) {
      ptr -> tse_id = ID2TSEIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REQUEST_REGET_TSE_details) {
      ptr -> details = ID2GetTSEDetailsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REQUEST_REGET_TSE_offset) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> offset = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2RequestReGetTSEFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2RequestReGetTSEAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2RequestReGetTSEAsnWrite(ID2RequestReGetTSEPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REQUEST_REGET_TSE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tse_id != NULL) {
      if ( ! ID2TSEIdAsnWrite(ptr -> tse_id, aip, ID2_REQUEST_REGET_TSE_tse_id)) {
         goto erret;
      }
   }
   if (ptr -> details != NULL) {
      if ( ! ID2GetTSEDetailsAsnWrite(ptr -> details, aip, ID2_REQUEST_REGET_TSE_details)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> offset;
   retval = AsnWrite(aip, ID2_REQUEST_REGET_TSE_offset,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SRequestGetChunksNew()
*
**************************************************/
NLM_EXTERN 
ID2SRequestGetChunksPtr LIBCALL
ID2SRequestGetChunksNew(void)
{
   ID2SRequestGetChunksPtr ptr = MemNew((size_t) sizeof(ID2SRequestGetChunks));

   return ptr;

}


/**************************************************
*
*    ID2SRequestGetChunksFree()
*
**************************************************/
NLM_EXTERN 
ID2SRequestGetChunksPtr LIBCALL
ID2SRequestGetChunksFree(ID2SRequestGetChunksPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2TSEIdFree(ptr -> tse_id);
   AsnGenericBaseSeqOfFree(ptr -> chunks ,ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SRequestGetChunksAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SRequestGetChunksPtr LIBCALL
ID2SRequestGetChunksAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SRequestGetChunksPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SRequestGetChunks ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_REQUEST_GET_CHUNKS);
   } else {
      atp = AsnLinkType(orig, ID2S_REQUEST_GET_CHUNKS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SRequestGetChunksNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_REQUEST_GET_CHUNKS_tse_id) {
      ptr -> tse_id = ID2TSEIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_REQUEST_GET_CHUNKS_chunks) {
      ptr -> chunks = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> chunks == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SRequestGetChunksFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SRequestGetChunksAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SRequestGetChunksAsnWrite(ID2SRequestGetChunksPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_REQUEST_GET_CHUNKS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tse_id != NULL) {
      if ( ! ID2TSEIdAsnWrite(ptr -> tse_id, aip, ID2S_REQUEST_GET_CHUNKS_tse_id)) {
         goto erret;
      }
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> chunks ,ASNCODE_INTVAL_SLOT, aip, ID2S_REQUEST_GET_CHUNKS_chunks, REQUEST_GET_CHUNKS_chunks_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2TSEIdNew()
*
**************************************************/
NLM_EXTERN 
ID2TSEIdPtr LIBCALL
ID2TSEIdNew(void)
{
   ID2TSEIdPtr ptr = MemNew((size_t) sizeof(ID2TSEId));

   return ptr;

}


/**************************************************
*
*    ID2TSEIdFree()
*
**************************************************/
NLM_EXTERN 
ID2TSEIdPtr LIBCALL
ID2TSEIdFree(ID2TSEIdPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    ID2TSEIdAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2TSEIdPtr LIBCALL
ID2TSEIdAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2TSEIdPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2TSEId ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_TSE_ID);
   } else {
      atp = AsnLinkType(orig, ID2_TSE_ID);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2TSEIdNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_TSE_ID_sat) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> sat = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_TSE_ID_sat_key) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> sat_key = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2TSEIdFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2TSEIdAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2TSEIdAsnWrite(ID2TSEIdPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_TSE_ID);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> sat;
   retval = AsnWrite(aip, ID2_TSE_ID_sat,  &av);
   av.intvalue = ptr -> sat_key;
   retval = AsnWrite(aip, ID2_TSE_ID_sat_key,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2GetTSEDetailsNew()
*
**************************************************/
NLM_EXTERN 
ID2GetTSEDetailsPtr LIBCALL
ID2GetTSEDetailsNew(void)
{
   ID2GetTSEDetailsPtr ptr = MemNew((size_t) sizeof(ID2GetTSEDetails));

   ptr -> seq_class_level = 1;
   ptr -> descr_level = 1;
   ptr -> descr_type_mask = 0;
   ptr -> annot_type_mask = 0;
   ptr -> feat_type_mask = 0;
   ptr -> sequence_level = 0;
   return ptr;

}


/**************************************************
*
*    ID2GetTSEDetailsFree()
*
**************************************************/
NLM_EXTERN 
ID2GetTSEDetailsPtr LIBCALL
ID2GetTSEDetailsFree(ID2GetTSEDetailsPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2SeqLocFree(ptr -> location);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2GetTSEDetailsAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2GetTSEDetailsPtr LIBCALL
ID2GetTSEDetailsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2GetTSEDetailsPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2GetTSEDetails ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_GET_TSE_DETAILS);
   } else {
      atp = AsnLinkType(orig, ID2_GET_TSE_DETAILS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2GetTSEDetailsNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_GET_TSE_DETAILS_location) {
      ptr -> location = ID2SeqLocAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TSE_DETAILS_seq_class_level) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> seq_class_level = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_GET_TSE_DETAILS_descr_level) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> descr_level = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TSE_DETAILS_descr_type_mask) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> descr_type_mask = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TSE_DETAILS_annot_type_mask) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> annot_type_mask = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GET_TSE_DETAILS_feat_type_mask) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> feat_type_mask = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GET_TSE_DETAILS_sequence_level) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> sequence_level = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2GetTSEDetailsFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2GetTSEDetailsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2GetTSEDetailsAsnWrite(ID2GetTSEDetailsPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_GET_TSE_DETAILS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> location != NULL) {
      if ( ! ID2SeqLocAsnWrite(ptr -> location, aip, ID2_GET_TSE_DETAILS_location)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> seq_class_level;
   retval = AsnWrite(aip, TSE_DETAILS_seq_class_level,  &av);
   av.intvalue = ptr -> descr_level;
   retval = AsnWrite(aip, ID2_GET_TSE_DETAILS_descr_level,  &av);
   av.intvalue = ptr -> descr_type_mask;
   retval = AsnWrite(aip, TSE_DETAILS_descr_type_mask,  &av);
   av.intvalue = ptr -> annot_type_mask;
   retval = AsnWrite(aip, TSE_DETAILS_annot_type_mask,  &av);
   av.intvalue = ptr -> feat_type_mask;
   retval = AsnWrite(aip, GET_TSE_DETAILS_feat_type_mask,  &av);
   av.intvalue = ptr -> sequence_level;
   retval = AsnWrite(aip, GET_TSE_DETAILS_sequence_level,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SeqLocFree()
*
**************************************************/
NLM_EXTERN 
ID2SeqLocPtr LIBCALL
ID2SeqLocFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ID2SeqLoc_int__:
      ID2IntervalFree(anp -> data.ptrvalue);
      break;
   case ID2SeqLoc_int_set:
      ID2PackedSeqIntsFree(anp -> data.ptrvalue);
      break;
   case ID2SeqLoc_whole_range:
      ID2IdRangeFree(anp -> data.ptrvalue);
      break;
   case ID2SeqLoc_loc_set:
      AsnGenericChoiceSeqOfFree((Pointer) pnt, (AsnOptFreeFunc) ID2SeqLocFree);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ID2SeqLocAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SeqLocPtr LIBCALL
ID2SeqLocAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SeqLoc ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_SEQ_LOC);
   } else {
      atp = AsnLinkType(orig, ID2_SEQ_LOC);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == ID2_SEQ_LOC_whole) {
      choice = ID2SeqLoc_whole;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == ID2_SEQ_LOC_int__) {
      choice = ID2SeqLoc_int__;
      func = (AsnReadFunc) ID2IntervalAsnRead;
   }
   else if (atp == ID2_SEQ_LOC_int_set) {
      choice = ID2SeqLoc_int_set;
      func = (AsnReadFunc) ID2PackedSeqIntsAsnRead;
   }
   else if (atp == ID2_SEQ_LOC_whole_range) {
      choice = ID2SeqLoc_whole_range;
      func = (AsnReadFunc) ID2IdRangeAsnRead;
   }
   else if (atp == ID2_SEQ_LOC_loc_set) {
      choice = ID2SeqLoc_loc_set;
      anp -> data.ptrvalue =
      AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2SeqLocAsnRead,             (AsnOptFreeFunc) ID2SeqLocFree);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ID2SeqLocAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SeqLocAsnWrite(ID2SeqLocPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ID2_SEQ_LOC);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ID2SeqLoc_whole:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, ID2_SEQ_LOC_whole, &av);
      break;
   case ID2SeqLoc_int__:
      writetype = ID2_SEQ_LOC_int__;
      func = (AsnWriteFunc) ID2IntervalAsnWrite;
      break;
   case ID2SeqLoc_int_set:
      writetype = ID2_SEQ_LOC_int_set;
      func = (AsnWriteFunc) ID2PackedSeqIntsAsnWrite;
      break;
   case ID2SeqLoc_whole_range:
      writetype = ID2_SEQ_LOC_whole_range;
      func = (AsnWriteFunc) ID2IdRangeAsnWrite;
      break;
   case ID2SeqLoc_loc_set:
      retval = AsnGenericChoiceSeqOfAsnWrite((Pointer) pnt, (AsnWriteFunc) ID2SeqLocAsnWrite, aip, ID2_SEQ_LOC_loc_set, ID2_SEQ_LOC_loc_set_E);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ID2ReplyNew()
*
**************************************************/
NLM_EXTERN 
ID2ReplyPtr LIBCALL
ID2ReplyNew(void)
{
   ID2ReplyPtr ptr = MemNew((size_t) sizeof(ID2Reply));

   return ptr;

}


/**************************************************
*
*    ID2ReplyFree()
*
**************************************************/
NLM_EXTERN 
ID2ReplyPtr LIBCALL
ID2ReplyFree(ID2ReplyPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2ParamsFree(ptr -> params);
   Reply_replyFree(ptr -> Reply_reply);
   AsnGenericUserSeqOfFree(ptr -> error, (AsnOptFreeFunc) ID2ErrorFree);
   return MemFree(ptr);
}


/**************************************************
*
*    Reply_replyFree()
*
**************************************************/
static 
Reply_replyPtr LIBCALL
Reply_replyFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case Reply_reply_get_package:
      ID2ReplyGetPackageFree(anp -> data.ptrvalue);
      break;
   case Reply_reply_seq_id_to_gi:
      ID2ReplySeqIdToGiFree(anp -> data.ptrvalue);
      break;
   case Reply_reply_gi_to_tse_id:
      ID2ReplyGiToTSEIdFree(anp -> data.ptrvalue);
      break;
   case Reply_reply_get_tse:
      ID2ReplyGetTSEFree(anp -> data.ptrvalue);
      break;
   case Reply_reply_get_tse_info:
      ID2SReplyGetTSEInfoFree(anp -> data.ptrvalue);
      break;
   case Reply_reply_get_chunk:
      ID2SReplyGetChunkFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ID2ReplyAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ReplyPtr LIBCALL
ID2ReplyAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ReplyPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2Reply ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REPLY);
   } else {
      atp = AsnLinkType(orig, ID2_REPLY);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ReplyNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REPLY_serial_number) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> serial_number = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_params) {
      ptr -> params = ID2ParamsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_reply) {
      ptr -> Reply_reply = Reply_replyAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_error) {
      ptr -> error = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2ErrorAsnRead, (AsnOptFreeFunc) ID2ErrorFree);
      if (isError && ptr -> error == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_end_of_reply) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> end_of_reply = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ReplyFree(ptr);
   goto ret;
}



/**************************************************
*
*    Reply_replyAsnRead()
*
**************************************************/
static 
Reply_replyPtr LIBCALL
Reply_replyAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Reply_reply ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REPLY_reply);
   } else {
      atp = AsnLinkType(orig, ID2_REPLY_reply);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == ID2_REPLY_reply_init) {
      choice = Reply_reply_init;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == ID2_REPLY_reply_get_package) {
      choice = Reply_reply_get_package;
      func = (AsnReadFunc) ID2ReplyGetPackageAsnRead;
   }
   else if (atp == ID2_REPLY_reply_seq_id_to_gi) {
      choice = Reply_reply_seq_id_to_gi;
      func = (AsnReadFunc) ID2ReplySeqIdToGiAsnRead;
   }
   else if (atp == ID2_REPLY_reply_gi_to_tse_id) {
      choice = Reply_reply_gi_to_tse_id;
      func = (AsnReadFunc) ID2ReplyGiToTSEIdAsnRead;
   }
   else if (atp == ID2_REPLY_reply_get_tse) {
      choice = Reply_reply_get_tse;
      func = (AsnReadFunc) ID2ReplyGetTSEAsnRead;
   }
   else if (atp == ID2_REPLY_reply_get_tse_info) {
      choice = Reply_reply_get_tse_info;
      func = (AsnReadFunc) ID2SReplyGetTSEInfoAsnRead;
   }
   else if (atp == ID2_REPLY_reply_get_chunk) {
      choice = Reply_reply_get_chunk;
      func = (AsnReadFunc) ID2SReplyGetChunkAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ID2ReplyAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ReplyAsnWrite(ID2ReplyPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REPLY);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> serial_number;
   retval = AsnWrite(aip, ID2_REPLY_serial_number,  &av);
   if (ptr -> params != NULL) {
      if ( ! ID2ParamsAsnWrite(ptr -> params, aip, ID2_REPLY_params)) {
         goto erret;
      }
   }
   if (ptr -> Reply_reply != NULL) {
      if ( ! Reply_replyAsnWrite(ptr -> Reply_reply, aip, ID2_REPLY_reply)) {
         goto erret;
      }
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> error, (AsnWriteFunc) ID2ErrorAsnWrite, aip, ID2_REPLY_error, ID2_REPLY_error_E);
   av.boolvalue = ptr -> end_of_reply;
   retval = AsnWrite(aip, ID2_REPLY_end_of_reply,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Reply_replyAsnWrite()
*
**************************************************/
static Boolean LIBCALL 
Reply_replyAsnWrite(Reply_replyPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ID2_REPLY_reply);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case Reply_reply_init:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, ID2_REPLY_reply_init, &av);
      break;
   case Reply_reply_get_package:
      writetype = ID2_REPLY_reply_get_package;
      func = (AsnWriteFunc) ID2ReplyGetPackageAsnWrite;
      break;
   case Reply_reply_seq_id_to_gi:
      writetype = ID2_REPLY_reply_seq_id_to_gi;
      func = (AsnWriteFunc) ID2ReplySeqIdToGiAsnWrite;
      break;
   case Reply_reply_gi_to_tse_id:
      writetype = ID2_REPLY_reply_gi_to_tse_id;
      func = (AsnWriteFunc) ID2ReplyGiToTSEIdAsnWrite;
      break;
   case Reply_reply_get_tse:
      writetype = ID2_REPLY_reply_get_tse;
      func = (AsnWriteFunc) ID2ReplyGetTSEAsnWrite;
      break;
   case Reply_reply_get_tse_info:
      writetype = ID2_REPLY_reply_get_tse_info;
      func = (AsnWriteFunc) ID2SReplyGetTSEInfoAsnWrite;
      break;
   case Reply_reply_get_chunk:
      writetype = ID2_REPLY_reply_get_chunk;
      func = (AsnWriteFunc) ID2SReplyGetChunkAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ID2ReplyGetPackageNew()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGetPackagePtr LIBCALL
ID2ReplyGetPackageNew(void)
{
   ID2ReplyGetPackagePtr ptr = MemNew((size_t) sizeof(ID2ReplyGetPackage));

   return ptr;

}


/**************************************************
*
*    ID2ReplyGetPackageFree()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGetPackagePtr LIBCALL
ID2ReplyGetPackageFree(ID2ReplyGetPackagePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   ID2ParamsFree(ptr -> params);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2ReplyGetPackageAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGetPackagePtr LIBCALL
ID2ReplyGetPackageAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ReplyGetPackagePtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2ReplyGetPackage ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REPLY_GET_PACKAGE);
   } else {
      atp = AsnLinkType(orig, ID2_REPLY_GET_PACKAGE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ReplyGetPackageNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REPLY_GET_PACKAGE_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_GET_PACKAGE_params) {
      ptr -> params = ID2ParamsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ReplyGetPackageFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ReplyGetPackageAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ReplyGetPackageAsnWrite(ID2ReplyGetPackagePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REPLY_GET_PACKAGE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, ID2_REPLY_GET_PACKAGE_name,  &av);
   }
   if (ptr -> params != NULL) {
      if ( ! ID2ParamsAsnWrite(ptr -> params, aip, ID2_REPLY_GET_PACKAGE_params)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2ReplySeqIdToGiNew()
*
**************************************************/
NLM_EXTERN 
ID2ReplySeqIdToGiPtr LIBCALL
ID2ReplySeqIdToGiNew(void)
{
   ID2ReplySeqIdToGiPtr ptr = MemNew((size_t) sizeof(ID2ReplySeqIdToGi));

   return ptr;

}


/**************************************************
*
*    ID2ReplySeqIdToGiFree()
*
**************************************************/
NLM_EXTERN 
ID2ReplySeqIdToGiPtr LIBCALL
ID2ReplySeqIdToGiFree(ID2ReplySeqIdToGiPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SeqIdFree(ptr -> seq_id);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2ReplySeqIdToGiAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ReplySeqIdToGiPtr LIBCALL
ID2ReplySeqIdToGiAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ReplySeqIdToGiPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2ReplySeqIdToGi ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REPLY_SEQ_ID_TO_GI);
   } else {
      atp = AsnLinkType(orig, ID2_REPLY_SEQ_ID_TO_GI);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ReplySeqIdToGiNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REPLY_SEQ_ID_TO_GI_seq_id) {
      ptr -> seq_id = SeqIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_SEQ_ID_TO_GI_gi) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gi = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ReplySeqIdToGiFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ReplySeqIdToGiAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ReplySeqIdToGiAsnWrite(ID2ReplySeqIdToGiPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REPLY_SEQ_ID_TO_GI);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> seq_id != NULL) {
      if ( ! SeqIdAsnWrite(ptr -> seq_id, aip, ID2_REPLY_SEQ_ID_TO_GI_seq_id)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> gi;
   retval = AsnWrite(aip, ID2_REPLY_SEQ_ID_TO_GI_gi,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2ReplyGiToTSEIdNew()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGiToTSEIdPtr LIBCALL
ID2ReplyGiToTSEIdNew(void)
{
   ID2ReplyGiToTSEIdPtr ptr = MemNew((size_t) sizeof(ID2ReplyGiToTSEId));

   ptr -> source = "0";
   return ptr;

}


/**************************************************
*
*    ID2ReplyGiToTSEIdFree()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGiToTSEIdPtr LIBCALL
ID2ReplyGiToTSEIdFree(ID2ReplyGiToTSEIdPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> source);
   AsnGenericUserSeqOfFree(ptr -> tses, (AsnOptFreeFunc) ID2TSEIdInfoFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2ReplyGiToTSEIdAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGiToTSEIdPtr LIBCALL
ID2ReplyGiToTSEIdAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ReplyGiToTSEIdPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2ReplyGiToTSEId ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REPLY_GI_TO_TSE_ID);
   } else {
      atp = AsnLinkType(orig, ID2_REPLY_GI_TO_TSE_ID);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ReplyGiToTSEIdNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REPLY_GI_TO_TSE_ID_gi) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gi = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_GI_TO_TSE_ID_source) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> source = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_GI_TO_TSE_ID_tses) {
      ptr -> tses = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2TSEIdInfoAsnRead, (AsnOptFreeFunc) ID2TSEIdInfoFree);
      if (isError && ptr -> tses == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ReplyGiToTSEIdFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ReplyGiToTSEIdAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ReplyGiToTSEIdAsnWrite(ID2ReplyGiToTSEIdPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REPLY_GI_TO_TSE_ID);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> gi;
   retval = AsnWrite(aip, ID2_REPLY_GI_TO_TSE_ID_gi,  &av);
   if (ptr -> source != NULL) {
      av.ptrvalue = ptr -> source;
      retval = AsnWrite(aip, ID2_REPLY_GI_TO_TSE_ID_source,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> tses, (AsnWriteFunc) ID2TSEIdInfoAsnWrite, aip, ID2_REPLY_GI_TO_TSE_ID_tses, ID2_REPLY_GI_TO_TSE_ID_tses_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2ReplyGetTSENew()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGetTSEPtr LIBCALL
ID2ReplyGetTSENew(void)
{
   ID2ReplyGetTSEPtr ptr = MemNew((size_t) sizeof(ID2ReplyGetTSE));

   return ptr;

}


/**************************************************
*
*    ID2ReplyGetTSEFree()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGetTSEPtr LIBCALL
ID2ReplyGetTSEFree(ID2ReplyGetTSEPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2TSEIdFree(ptr -> tse_id);
   ID2ReplyDataFree(ptr -> data);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2ReplyGetTSEAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ReplyGetTSEPtr LIBCALL
ID2ReplyGetTSEAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ReplyGetTSEPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2ReplyGetTSE ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REPLY_GET_TSE);
   } else {
      atp = AsnLinkType(orig, ID2_REPLY_GET_TSE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ReplyGetTSENew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REPLY_GET_TSE_tse_id) {
      ptr -> tse_id = ID2TSEIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_GET_TSE_data) {
      ptr -> data = ID2ReplyDataAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ReplyGetTSEFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ReplyGetTSEAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ReplyGetTSEAsnWrite(ID2ReplyGetTSEPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REPLY_GET_TSE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tse_id != NULL) {
      if ( ! ID2TSEIdAsnWrite(ptr -> tse_id, aip, ID2_REPLY_GET_TSE_tse_id)) {
         goto erret;
      }
   }
   if (ptr -> data != NULL) {
      if ( ! ID2ReplyDataAsnWrite(ptr -> data, aip, ID2_REPLY_GET_TSE_data)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SReplyGetTSEInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SReplyGetTSEInfoPtr LIBCALL
ID2SReplyGetTSEInfoNew(void)
{
   ID2SReplyGetTSEInfoPtr ptr = MemNew((size_t) sizeof(ID2SReplyGetTSEInfo));

   return ptr;

}


/**************************************************
*
*    ID2SReplyGetTSEInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SReplyGetTSEInfoPtr LIBCALL
ID2SReplyGetTSEInfoFree(ID2SReplyGetTSEInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2TSEIdFree(ptr -> tse_id);
   ID2ReplyDataFree(ptr -> info);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SReplyGetTSEInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SReplyGetTSEInfoPtr LIBCALL
ID2SReplyGetTSEInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SReplyGetTSEInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SReplyGetTSEInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_REPLY_GET_TSE_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_REPLY_GET_TSE_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SReplyGetTSEInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_REPLY_GET_TSE_INFO_tse_id) {
      ptr -> tse_id = ID2TSEIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GET_TSE_INFO_split_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> split_version = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_REPLY_GET_TSE_INFO_info) {
      ptr -> info = ID2ReplyDataAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SReplyGetTSEInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SReplyGetTSEInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SReplyGetTSEInfoAsnWrite(ID2SReplyGetTSEInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_REPLY_GET_TSE_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tse_id != NULL) {
      if ( ! ID2TSEIdAsnWrite(ptr -> tse_id, aip, ID2S_REPLY_GET_TSE_INFO_tse_id)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> split_version;
   retval = AsnWrite(aip, GET_TSE_INFO_split_version,  &av);
   if (ptr -> info != NULL) {
      if ( ! ID2ReplyDataAsnWrite(ptr -> info, aip, ID2S_REPLY_GET_TSE_INFO_info)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SReplyGetChunkNew()
*
**************************************************/
NLM_EXTERN 
ID2SReplyGetChunkPtr LIBCALL
ID2SReplyGetChunkNew(void)
{
   ID2SReplyGetChunkPtr ptr = MemNew((size_t) sizeof(ID2SReplyGetChunk));

   return ptr;

}


/**************************************************
*
*    ID2SReplyGetChunkFree()
*
**************************************************/
NLM_EXTERN 
ID2SReplyGetChunkPtr LIBCALL
ID2SReplyGetChunkFree(ID2SReplyGetChunkPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2TSEIdFree(ptr -> tse_id);
   ID2ReplyDataFree(ptr -> data);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SReplyGetChunkAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SReplyGetChunkPtr LIBCALL
ID2SReplyGetChunkAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SReplyGetChunkPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SReplyGetChunk ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_REPLY_GET_CHUNK);
   } else {
      atp = AsnLinkType(orig, ID2S_REPLY_GET_CHUNK);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SReplyGetChunkNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_REPLY_GET_CHUNK_tse_id) {
      ptr -> tse_id = ID2TSEIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_REPLY_GET_CHUNK_chunk_id) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> chunk_id = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_REPLY_GET_CHUNK_data) {
      ptr -> data = ID2ReplyDataAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SReplyGetChunkFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SReplyGetChunkAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SReplyGetChunkAsnWrite(ID2SReplyGetChunkPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_REPLY_GET_CHUNK);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tse_id != NULL) {
      if ( ! ID2TSEIdAsnWrite(ptr -> tse_id, aip, ID2S_REPLY_GET_CHUNK_tse_id)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> chunk_id;
   retval = AsnWrite(aip, ID2S_REPLY_GET_CHUNK_chunk_id,  &av);
   if (ptr -> data != NULL) {
      if ( ! ID2ReplyDataAsnWrite(ptr -> data, aip, ID2S_REPLY_GET_CHUNK_data)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2ErrorNew()
*
**************************************************/
NLM_EXTERN 
ID2ErrorPtr LIBCALL
ID2ErrorNew(void)
{
   ID2ErrorPtr ptr = MemNew((size_t) sizeof(ID2Error));

   return ptr;

}


/**************************************************
*
*    ID2ErrorFree()
*
**************************************************/
NLM_EXTERN 
ID2ErrorPtr LIBCALL
ID2ErrorFree(ID2ErrorPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> message);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2ErrorAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ErrorPtr LIBCALL
ID2ErrorAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ErrorPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2Error ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_ERROR);
   } else {
      atp = AsnLinkType(orig, ID2_ERROR);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ErrorNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_ERROR_severity) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> severity = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_ERROR_retry_delay) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> retry_delay = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_ERROR_message) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> message = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ErrorFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ErrorAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ErrorAsnWrite(ID2ErrorPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_ERROR);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> severity;
   retval = AsnWrite(aip, ID2_ERROR_severity,  &av);
   av.intvalue = ptr -> retry_delay;
   retval = AsnWrite(aip, ID2_ERROR_retry_delay,  &av);
   if (ptr -> message != NULL) {
      av.ptrvalue = ptr -> message;
      retval = AsnWrite(aip, ID2_ERROR_message,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2TSEIdInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2TSEIdInfoPtr LIBCALL
ID2TSEIdInfoNew(void)
{
   ID2TSEIdInfoPtr ptr = MemNew((size_t) sizeof(ID2TSEIdInfo));

   ptr -> split_version = 0;
   return ptr;

}


/**************************************************
*
*    ID2TSEIdInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2TSEIdInfoPtr LIBCALL
ID2TSEIdInfoFree(ID2TSEIdInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2TSEIdFree(ptr -> tse_id);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2TSEIdInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2TSEIdInfoPtr LIBCALL
ID2TSEIdInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2TSEIdInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2TSEIdInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_TSE_ID_INFO);
   } else {
      atp = AsnLinkType(orig, ID2_TSE_ID_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2TSEIdInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_TSE_ID_INFO_tse_id) {
      ptr -> tse_id = ID2TSEIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_TSE_ID_INFO_split_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> split_version = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2TSEIdInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2TSEIdInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2TSEIdInfoAsnWrite(ID2TSEIdInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_TSE_ID_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tse_id != NULL) {
      if ( ! ID2TSEIdAsnWrite(ptr -> tse_id, aip, ID2_TSE_ID_INFO_tse_id)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> split_version;
   retval = AsnWrite(aip, ID2_TSE_ID_INFO_split_version,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2ReplyDataNew()
*
**************************************************/
NLM_EXTERN 
ID2ReplyDataPtr LIBCALL
ID2ReplyDataNew(void)
{
   ID2ReplyDataPtr ptr = MemNew((size_t) sizeof(ID2ReplyData));

   return ptr;

}


/**************************************************
*
*    ID2ReplyDataFree()
*
**************************************************/
NLM_EXTERN 
ID2ReplyDataPtr LIBCALL
ID2ReplyDataFree(ID2ReplyDataPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> data ,ASNCODE_BYTEVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2ReplyDataAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ReplyDataPtr LIBCALL
ID2ReplyDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ReplyDataPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2ReplyData ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_REPLY_DATA);
   } else {
      atp = AsnLinkType(orig, ID2_REPLY_DATA);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ReplyDataNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_REPLY_DATA_data_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> data_type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_DATA_data_format) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> data_format = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_DATA_data_compression) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> data_compression = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_REPLY_DATA_data) {
      ptr -> data = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_BYTEVAL_SLOT, &isError);
      if (isError && ptr -> data == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ReplyDataFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ReplyDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ReplyDataAsnWrite(ID2ReplyDataPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_REPLY_DATA);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> data_type;
   retval = AsnWrite(aip, ID2_REPLY_DATA_data_type,  &av);
   av.intvalue = ptr -> data_format;
   retval = AsnWrite(aip, ID2_REPLY_DATA_data_format,  &av);
   av.intvalue = ptr -> data_compression;
   retval = AsnWrite(aip, ID2_REPLY_DATA_data_compression,  &av);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> data ,ASNCODE_BYTEVAL_SLOT, aip, ID2_REPLY_DATA_data, ID2_REPLY_DATA_data_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SSplitInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SSplitInfoPtr LIBCALL
ID2SSplitInfoNew(void)
{
   ID2SSplitInfoPtr ptr = MemNew((size_t) sizeof(ID2SSplitInfo));

   return ptr;

}


/**************************************************
*
*    ID2SSplitInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SSplitInfoPtr LIBCALL
ID2SSplitInfoFree(ID2SSplitInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr -> bioseqs_info, (AsnOptFreeFunc) ID2SBioseqsInfoFree);
   AsnGenericUserSeqOfFree(ptr -> chunks, (AsnOptFreeFunc) ID2SChunkInfoFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SSplitInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SSplitInfoPtr LIBCALL
ID2SSplitInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SSplitInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SSplitInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_SPLIT_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_SPLIT_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SSplitInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_SPLIT_INFO_bioseqs_info) {
      ptr -> bioseqs_info = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2SBioseqsInfoAsnRead, (AsnOptFreeFunc) ID2SBioseqsInfoFree);
      if (isError && ptr -> bioseqs_info == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_SPLIT_INFO_chunks) {
      ptr -> chunks = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2SChunkInfoAsnRead, (AsnOptFreeFunc) ID2SChunkInfoFree);
      if (isError && ptr -> chunks == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SSplitInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SSplitInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SSplitInfoAsnWrite(ID2SSplitInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_SPLIT_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   AsnGenericUserSeqOfAsnWrite(ptr -> bioseqs_info, (AsnWriteFunc) ID2SBioseqsInfoAsnWrite, aip, ID2S_SPLIT_INFO_bioseqs_info, ID2S_SPLIT_INFO_bioseqs_info_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> chunks, (AsnWriteFunc) ID2SChunkInfoAsnWrite, aip, ID2S_SPLIT_INFO_chunks, ID2S_SPLIT_INFO_chunks_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SBioseqsInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SBioseqsInfoPtr LIBCALL
ID2SBioseqsInfoNew(void)
{
   ID2SBioseqsInfoPtr ptr = MemNew((size_t) sizeof(ID2SBioseqsInfo));

   return ptr;

}


/**************************************************
*
*    ID2SBioseqsInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SBioseqsInfoPtr LIBCALL
ID2SBioseqsInfoFree(ID2SBioseqsInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2SBioseqInfoFree(ptr -> info);
   ID2IdRangeFree(ptr -> bioseqs);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SBioseqsInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SBioseqsInfoPtr LIBCALL
ID2SBioseqsInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SBioseqsInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SBioseqsInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_BIOSEQS_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_BIOSEQS_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SBioseqsInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_BIOSEQS_INFO_info) {
      ptr -> info = ID2SBioseqInfoAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_BIOSEQS_INFO_bioseqs) {
      ptr -> bioseqs = ID2IdRangeAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SBioseqsInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SBioseqsInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SBioseqsInfoAsnWrite(ID2SBioseqsInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_BIOSEQS_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> info != NULL) {
      if ( ! ID2SBioseqInfoAsnWrite(ptr -> info, aip, ID2S_BIOSEQS_INFO_info)) {
         goto erret;
      }
   }
   if (ptr -> bioseqs != NULL) {
      if ( ! ID2IdRangeAsnWrite(ptr -> bioseqs, aip, ID2S_BIOSEQS_INFO_bioseqs)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SChunkInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SChunkInfoPtr LIBCALL
ID2SChunkInfoNew(void)
{
   ID2SChunkInfoPtr ptr = MemNew((size_t) sizeof(ID2SChunkInfo));

   return ptr;

}


/**************************************************
*
*    ID2SChunkInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SChunkInfoPtr LIBCALL
ID2SChunkInfoFree(ID2SChunkInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr -> content, (AsnOptFreeFunc) ID2SChunkContentFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SChunkInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SChunkInfoPtr LIBCALL
ID2SChunkInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SChunkInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SChunkInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_CHUNK_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_CHUNK_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SChunkInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_CHUNK_INFO_id) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> id = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_CHUNK_INFO_content) {
      ptr -> content = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2SChunkContentAsnRead, (AsnOptFreeFunc) ID2SChunkContentFree);
      if (isError && ptr -> content == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SChunkInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SChunkInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SChunkInfoAsnWrite(ID2SChunkInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_CHUNK_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> id;
   retval = AsnWrite(aip, ID2S_CHUNK_INFO_id,  &av);
   AsnGenericChoiceSeqOfAsnWrite(ptr -> content, (AsnWriteFunc) ID2SChunkContentAsnWrite, aip, ID2S_CHUNK_INFO_content, ID2S_CHUNK_INFO_content_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SBioseqInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SBioseqInfoPtr LIBCALL
ID2SBioseqInfoNew(void)
{
   ID2SBioseqInfoPtr ptr = MemNew((size_t) sizeof(ID2SBioseqInfo));

   return ptr;

}


/**************************************************
*
*    ID2SBioseqInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SBioseqInfoPtr LIBCALL
ID2SBioseqInfoFree(ID2SBioseqInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ID2SSequenceSplitInfoFree(ptr -> sequence_split);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SBioseqInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SBioseqInfoPtr LIBCALL
ID2SBioseqInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SBioseqInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SBioseqInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_BIOSEQ_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_BIOSEQ_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SBioseqInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_BIOSEQ_INFO_gap_count) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gap_count = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSEQ_INFO_seq_map_has_ref) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> seq_map_has_ref = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_BIOSEQ_INFO_sequence_split) {
      ptr -> sequence_split = ID2SSequenceSplitInfoAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SBioseqInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SBioseqInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SBioseqInfoAsnWrite(ID2SBioseqInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_BIOSEQ_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> gap_count;
   retval = AsnWrite(aip, ID2S_BIOSEQ_INFO_gap_count,  &av);
   av.boolvalue = ptr -> seq_map_has_ref;
   retval = AsnWrite(aip, BIOSEQ_INFO_seq_map_has_ref,  &av);
   if (ptr -> sequence_split != NULL) {
      if ( ! ID2SSequenceSplitInfoAsnWrite(ptr -> sequence_split, aip, ID2S_BIOSEQ_INFO_sequence_split)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2IdRangeNew()
*
**************************************************/
NLM_EXTERN 
ID2IdRangePtr LIBCALL
ID2IdRangeNew(void)
{
   ID2IdRangePtr ptr = MemNew((size_t) sizeof(ID2IdRange));

   ptr -> count = 1;
   return ptr;

}


/**************************************************
*
*    ID2IdRangeFree()
*
**************************************************/
NLM_EXTERN 
ID2IdRangePtr LIBCALL
ID2IdRangeFree(ID2IdRangePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    ID2IdRangeAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2IdRangePtr LIBCALL
ID2IdRangeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2IdRangePtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2IdRange ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_ID_RANGE);
   } else {
      atp = AsnLinkType(orig, ID2_ID_RANGE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2IdRangeNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_ID_RANGE_start) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> start = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_ID_RANGE_count) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> count = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2IdRangeFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2IdRangeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2IdRangeAsnWrite(ID2IdRangePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_ID_RANGE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> start;
   retval = AsnWrite(aip, ID2_ID_RANGE_start,  &av);
   av.intvalue = ptr -> count;
   retval = AsnWrite(aip, ID2_ID_RANGE_count,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SSequenceSplitInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SSequenceSplitInfoPtr LIBCALL
ID2SSequenceSplitInfoNew(void)
{
   ID2SSequenceSplitInfoPtr ptr = MemNew((size_t) sizeof(ID2SSequenceSplitInfo));

   return ptr;

}


/**************************************************
*
*    ID2SSequenceSplitInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SSequenceSplitInfoPtr LIBCALL
ID2SSequenceSplitInfoFree(ID2SSequenceSplitInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> chunk_blocks ,ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SSequenceSplitInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SSequenceSplitInfoPtr LIBCALL
ID2SSequenceSplitInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SSequenceSplitInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SSequenceSplitInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_SEQUENCE_SPLIT_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_SEQUENCE_SPLIT_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SSequenceSplitInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SEQUENCE_SPLIT_INFO_block_size) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> block_size = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SPLIT_INFO_chunk_start) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> chunk_start = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SPLIT_INFO_chunk_blocks) {
      ptr -> chunk_blocks = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> chunk_blocks == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SSequenceSplitInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SSequenceSplitInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SSequenceSplitInfoAsnWrite(ID2SSequenceSplitInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_SEQUENCE_SPLIT_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> block_size;
   retval = AsnWrite(aip, SEQUENCE_SPLIT_INFO_block_size,  &av);
   av.intvalue = ptr -> chunk_start;
   retval = AsnWrite(aip, SPLIT_INFO_chunk_start,  &av);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> chunk_blocks ,ASNCODE_INTVAL_SLOT, aip, SPLIT_INFO_chunk_blocks, SPLIT_INFO_chunk_blocks_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SChunkContentFree()
*
**************************************************/
NLM_EXTERN 
ID2SChunkContentPtr LIBCALL
ID2SChunkContentFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case ID2SChunkContent_seq_descr:
      ID2SSeqDescrInfoFree(anp -> data.ptrvalue);
      break;
   case ID2SChunkContent_seq_annot:
      ID2SSeqAnnotInfoFree(anp -> data.ptrvalue);
      break;
   case ID2SChunkContent_seq_assembly:
      ID2SSeqAssemblyInfoFree(anp -> data.ptrvalue);
      break;
   case ID2SChunkContent_seq_map:
      ID2SeqLocFree(anp -> data.ptrvalue);
      break;
   case ID2SChunkContent_seq_data:
      ID2SeqLocFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ID2SChunkContentAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SChunkContentPtr LIBCALL
ID2SChunkContentAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SChunkContent ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_CHUNK_CONTENT);
   } else {
      atp = AsnLinkType(orig, ID2S_CHUNK_CONTENT);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == ID2S_CHUNK_CONTENT_seq_descr) {
      choice = ID2SChunkContent_seq_descr;
      func = (AsnReadFunc) ID2SSeqDescrInfoAsnRead;
   }
   else if (atp == ID2S_CHUNK_CONTENT_seq_annot) {
      choice = ID2SChunkContent_seq_annot;
      func = (AsnReadFunc) ID2SSeqAnnotInfoAsnRead;
   }
   else if (atp == ID2S_CHUNK_CONTENT_seq_assembly) {
      choice = ID2SChunkContent_seq_assembly;
      func = (AsnReadFunc) ID2SSeqAssemblyInfoAsnRead;
   }
   else if (atp == ID2S_CHUNK_CONTENT_seq_map) {
      choice = ID2SChunkContent_seq_map;
      func = (AsnReadFunc) ID2SeqLocAsnRead;
   }
   else if (atp == ID2S_CHUNK_CONTENT_seq_data) {
      choice = ID2SChunkContent_seq_data;
      func = (AsnReadFunc) ID2SeqLocAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ID2SChunkContentAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SChunkContentAsnWrite(ID2SChunkContentPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ID2S_CHUNK_CONTENT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case ID2SChunkContent_seq_descr:
      writetype = ID2S_CHUNK_CONTENT_seq_descr;
      func = (AsnWriteFunc) ID2SSeqDescrInfoAsnWrite;
      break;
   case ID2SChunkContent_seq_annot:
      writetype = ID2S_CHUNK_CONTENT_seq_annot;
      func = (AsnWriteFunc) ID2SSeqAnnotInfoAsnWrite;
      break;
   case ID2SChunkContent_seq_assembly:
      writetype = ID2S_CHUNK_CONTENT_seq_assembly;
      func = (AsnWriteFunc) ID2SSeqAssemblyInfoAsnWrite;
      break;
   case ID2SChunkContent_seq_map:
      writetype = ID2S_CHUNK_CONTENT_seq_map;
      func = (AsnWriteFunc) ID2SeqLocAsnWrite;
      break;
   case ID2SChunkContent_seq_data:
      writetype = ID2S_CHUNK_CONTENT_seq_data;
      func = (AsnWriteFunc) ID2SeqLocAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ID2SSeqDescrInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SSeqDescrInfoPtr LIBCALL
ID2SSeqDescrInfoNew(void)
{
   ID2SSeqDescrInfoPtr ptr = MemNew((size_t) sizeof(ID2SSeqDescrInfo));

   return ptr;

}


/**************************************************
*
*    ID2SSeqDescrInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SSeqDescrInfoPtr LIBCALL
ID2SSeqDescrInfoFree(ID2SSeqDescrInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr -> bioseqs, (AsnOptFreeFunc) ID2IdRangeFree);
   AsnGenericUserSeqOfFree(ptr -> bioseq_sets, (AsnOptFreeFunc) ID2IdRangeFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SSeqDescrInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SSeqDescrInfoPtr LIBCALL
ID2SSeqDescrInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SSeqDescrInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SSeqDescrInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_SEQ_DESCR_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_SEQ_DESCR_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SSeqDescrInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_SEQ_DESCR_INFO_type_mask) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type_mask = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_SEQ_DESCR_INFO_bioseqs) {
      ptr -> bioseqs = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2IdRangeAsnRead, (AsnOptFreeFunc) ID2IdRangeFree);
      if (isError && ptr -> bioseqs == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_SEQ_DESCR_INFO_bioseq_sets) {
      ptr -> bioseq_sets = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2IdRangeAsnRead, (AsnOptFreeFunc) ID2IdRangeFree);
      if (isError && ptr -> bioseq_sets == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SSeqDescrInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SSeqDescrInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SSeqDescrInfoAsnWrite(ID2SSeqDescrInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_SEQ_DESCR_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type_mask;
   retval = AsnWrite(aip, ID2S_SEQ_DESCR_INFO_type_mask,  &av);
   AsnGenericUserSeqOfAsnWrite(ptr -> bioseqs, (AsnWriteFunc) ID2IdRangeAsnWrite, aip, ID2S_SEQ_DESCR_INFO_bioseqs, ID2S_SEQ_DESCR_INFO_bioseqs_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> bioseq_sets, (AsnWriteFunc) ID2IdRangeAsnWrite, aip, ID2S_SEQ_DESCR_INFO_bioseq_sets, SEQ_DESCR_INFO_bioseq_sets_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SSeqAnnotInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SSeqAnnotInfoPtr LIBCALL
ID2SSeqAnnotInfoNew(void)
{
   ID2SSeqAnnotInfoPtr ptr = MemNew((size_t) sizeof(ID2SSeqAnnotInfo));

   return ptr;

}


/**************************************************
*
*    ID2SSeqAnnotInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SSeqAnnotInfoPtr LIBCALL
ID2SSeqAnnotInfoFree(ID2SSeqAnnotInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   AsnGenericUserSeqOfFree(ptr -> feat, (AsnOptFreeFunc) ID2SFeatTypeInfoFree);
   ID2SeqLocFree(ptr -> seq_loc);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SSeqAnnotInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SSeqAnnotInfoPtr LIBCALL
ID2SSeqAnnotInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SSeqAnnotInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SSeqAnnotInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_SEQ_ANNOT_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_SEQ_ANNOT_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SSeqAnnotInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_SEQ_ANNOT_INFO_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_SEQ_ANNOT_INFO_align) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> align = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_SEQ_ANNOT_INFO_graph) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> graph = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_SEQ_ANNOT_INFO_feat) {
      ptr -> feat = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2SFeatTypeInfoAsnRead, (AsnOptFreeFunc) ID2SFeatTypeInfoFree);
      if (isError && ptr -> feat == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_SEQ_ANNOT_INFO_seq_loc) {
      ptr -> seq_loc = ID2SeqLocAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SSeqAnnotInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SSeqAnnotInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SSeqAnnotInfoAsnWrite(ID2SSeqAnnotInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_SEQ_ANNOT_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, ID2S_SEQ_ANNOT_INFO_name,  &av);
   }
   av.boolvalue = ptr -> align;
   retval = AsnWrite(aip, ID2S_SEQ_ANNOT_INFO_align,  &av);
   av.boolvalue = ptr -> graph;
   retval = AsnWrite(aip, ID2S_SEQ_ANNOT_INFO_graph,  &av);
   AsnGenericUserSeqOfAsnWrite(ptr -> feat, (AsnWriteFunc) ID2SFeatTypeInfoAsnWrite, aip, ID2S_SEQ_ANNOT_INFO_feat, ID2S_SEQ_ANNOT_INFO_feat_E);
   if (ptr -> seq_loc != NULL) {
      if ( ! ID2SeqLocAsnWrite(ptr -> seq_loc, aip, ID2S_SEQ_ANNOT_INFO_seq_loc)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SSeqAssemblyInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SSeqAssemblyInfoPtr LIBCALL
ID2SSeqAssemblyInfoNew(void)
{
   ID2SSeqAssemblyInfoPtr ptr = MemNew((size_t) sizeof(ID2SSeqAssemblyInfo));

   return ptr;

}


/**************************************************
*
*    ID2SSeqAssemblyInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SSeqAssemblyInfoPtr LIBCALL
ID2SSeqAssemblyInfoFree(ID2SSeqAssemblyInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr -> bioseqs, (AsnOptFreeFunc) ID2IdRangeFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SSeqAssemblyInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SSeqAssemblyInfoPtr LIBCALL
ID2SSeqAssemblyInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SSeqAssemblyInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SSeqAssemblyInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_SEQ_ASSEMBLY_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_SEQ_ASSEMBLY_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SSeqAssemblyInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_SEQ_ASSEMBLY_INFO_bioseqs) {
      ptr -> bioseqs = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2IdRangeAsnRead, (AsnOptFreeFunc) ID2IdRangeFree);
      if (isError && ptr -> bioseqs == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SSeqAssemblyInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SSeqAssemblyInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SSeqAssemblyInfoAsnWrite(ID2SSeqAssemblyInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_SEQ_ASSEMBLY_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   AsnGenericUserSeqOfAsnWrite(ptr -> bioseqs, (AsnWriteFunc) ID2IdRangeAsnWrite, aip, ID2S_SEQ_ASSEMBLY_INFO_bioseqs, SEQ_ASSEMBLY_INFO_bioseqs_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SFeatTypeInfoNew()
*
**************************************************/
NLM_EXTERN 
ID2SFeatTypeInfoPtr LIBCALL
ID2SFeatTypeInfoNew(void)
{
   ID2SFeatTypeInfoPtr ptr = MemNew((size_t) sizeof(ID2SFeatTypeInfo));

   return ptr;

}


/**************************************************
*
*    ID2SFeatTypeInfoFree()
*
**************************************************/
NLM_EXTERN 
ID2SFeatTypeInfoPtr LIBCALL
ID2SFeatTypeInfoFree(ID2SFeatTypeInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> subtypes ,ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SFeatTypeInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SFeatTypeInfoPtr LIBCALL
ID2SFeatTypeInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SFeatTypeInfoPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SFeatTypeInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_FEAT_TYPE_INFO);
   } else {
      atp = AsnLinkType(orig, ID2S_FEAT_TYPE_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SFeatTypeInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_FEAT_TYPE_INFO_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_FEAT_TYPE_INFO_subtypes) {
      ptr -> subtypes = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> subtypes == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SFeatTypeInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SFeatTypeInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SFeatTypeInfoAsnWrite(ID2SFeatTypeInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_FEAT_TYPE_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, ID2S_FEAT_TYPE_INFO_type,  &av);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> subtypes ,ASNCODE_INTVAL_SLOT, aip, ID2S_FEAT_TYPE_INFO_subtypes, ID2S_FEAT_TYPE_INFO_subtypes_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SChunkNew()
*
**************************************************/
NLM_EXTERN 
ID2SChunkPtr LIBCALL
ID2SChunkNew(void)
{
   ID2SChunkPtr ptr = MemNew((size_t) sizeof(ID2SChunk));

   return ptr;

}


/**************************************************
*
*    ID2SChunkFree()
*
**************************************************/
NLM_EXTERN 
ID2SChunkPtr LIBCALL
ID2SChunkFree(ID2SChunkPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr -> data, (AsnOptFreeFunc) ID2SChunkDataFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SChunkAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SChunkPtr LIBCALL
ID2SChunkAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SChunkPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SChunk ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_CHUNK);
   } else {
      atp = AsnLinkType(orig, ID2S_CHUNK);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SChunkNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_CHUNK_data) {
      ptr -> data = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2SChunkDataAsnRead, (AsnOptFreeFunc) ID2SChunkDataFree);
      if (isError && ptr -> data == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SChunkFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SChunkAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SChunkAsnWrite(ID2SChunkPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_CHUNK);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   AsnGenericUserSeqOfAsnWrite(ptr -> data, (AsnWriteFunc) ID2SChunkDataAsnWrite, aip, ID2S_CHUNK_data, ID2S_CHUNK_data_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SChunkDataNew()
*
**************************************************/
NLM_EXTERN 
ID2SChunkDataPtr LIBCALL
ID2SChunkDataNew(void)
{
   ID2SChunkDataPtr ptr = MemNew((size_t) sizeof(ID2SChunkData));

   return ptr;

}


/**************************************************
*
*    ID2SChunkDataFree()
*
**************************************************/
NLM_EXTERN 
ID2SChunkDataPtr LIBCALL
ID2SChunkDataFree(ID2SChunkDataPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   Id_idFree(ptr -> Id_id);
   AsnGenericChoiceSeqOfFree(ptr -> descrs, (AsnOptFreeFunc) SeqDescrFree);
   AsnGenericUserSeqOfFree(ptr -> annots, (AsnOptFreeFunc) SeqAnnotFree);
   AsnGenericUserSeqOfFree(ptr -> assembly, (AsnOptFreeFunc) SeqAlignFree);
   AsnGenericUserSeqOfFree(ptr -> seq_map, (AsnOptFreeFunc) SeqLiteralFree);
   AsnGenericUserSeqOfFree(ptr -> seq_data, (AsnOptFreeFunc) SeqLiteralFree);
   return MemFree(ptr);
}


/**************************************************
*
*    Id_idFree()
*
**************************************************/
static 
Id_idPtr LIBCALL
Id_idFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    ID2SChunkDataAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SChunkDataPtr LIBCALL
ID2SChunkDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SChunkDataPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SChunkData ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_CHUNK_DATA);
   } else {
      atp = AsnLinkType(orig, ID2S_CHUNK_DATA);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SChunkDataNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2S_CHUNK_DATA_id) {
      ptr -> Id_id = Id_idAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_CHUNK_DATA_descrs) {
      ptr -> descrs = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqDescrAsnRead, (AsnOptFreeFunc) SeqDescrFree);
      if (isError && ptr -> descrs == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_CHUNK_DATA_annots) {
      ptr -> annots = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqAnnotAsnRead, (AsnOptFreeFunc) SeqAnnotFree);
      if (isError && ptr -> annots == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_CHUNK_DATA_assembly) {
      ptr -> assembly = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqAlignAsnRead, (AsnOptFreeFunc) SeqAlignFree);
      if (isError && ptr -> assembly == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_CHUNK_DATA_seq_map) {
      ptr -> seq_map = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqLiteralAsnRead, (AsnOptFreeFunc) SeqLiteralFree);
      if (isError && ptr -> seq_map == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2S_CHUNK_DATA_seq_data) {
      ptr -> seq_data = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqLiteralAsnRead, (AsnOptFreeFunc) SeqLiteralFree);
      if (isError && ptr -> seq_data == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SChunkDataFree(ptr);
   goto ret;
}



/**************************************************
*
*    Id_idAsnRead()
*
**************************************************/
static 
Id_idPtr LIBCALL
Id_idAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Id_id ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2S_CHUNK_DATA_id);
   } else {
      atp = AsnLinkType(orig, ID2S_CHUNK_DATA_id);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == ID2S_CHUNK_DATA_id_bioseq_set) {
      choice = Id_id_bioseq_set;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == ID2S_CHUNK_DATA_id_gi) {
      choice = Id_id_gi;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    ID2SChunkDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SChunkDataAsnWrite(ID2SChunkDataPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2S_CHUNK_DATA);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> Id_id != NULL) {
      if ( ! Id_idAsnWrite(ptr -> Id_id, aip, ID2S_CHUNK_DATA_id)) {
         goto erret;
      }
   }
   AsnGenericChoiceSeqOfAsnWrite(ptr -> descrs, (AsnWriteFunc) SeqDescrAsnWrite, aip, ID2S_CHUNK_DATA_descrs, ID2S_CHUNK_DATA_descrs_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> annots, (AsnWriteFunc) SeqAnnotAsnWrite, aip, ID2S_CHUNK_DATA_annots, ID2S_CHUNK_DATA_annots_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> assembly, (AsnWriteFunc) SeqAlignAsnWrite, aip, ID2S_CHUNK_DATA_assembly, ID2S_CHUNK_DATA_assembly_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> seq_map, (AsnWriteFunc) SeqLiteralAsnWrite, aip, ID2S_CHUNK_DATA_seq_map, ID2S_CHUNK_DATA_seq_map_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> seq_data, (AsnWriteFunc) SeqLiteralAsnWrite, aip, ID2S_CHUNK_DATA_seq_data, ID2S_CHUNK_DATA_seq_data_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Id_idAsnWrite()
*
**************************************************/
static Boolean LIBCALL 
Id_idAsnWrite(Id_idPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ID2S_CHUNK_DATA_id);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case Id_id_bioseq_set:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, ID2S_CHUNK_DATA_id_bioseq_set, &av);
      break;
   case Id_id_gi:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, ID2S_CHUNK_DATA_id_gi, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    ID2IntervalNew()
*
**************************************************/
NLM_EXTERN 
ID2IntervalPtr LIBCALL
ID2IntervalNew(void)
{
   ID2IntervalPtr ptr = MemNew((size_t) sizeof(ID2Interval));

   ptr -> length = 1;
   return ptr;

}


/**************************************************
*
*    ID2IntervalFree()
*
**************************************************/
NLM_EXTERN 
ID2IntervalPtr LIBCALL
ID2IntervalFree(ID2IntervalPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    ID2IntervalAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2IntervalPtr LIBCALL
ID2IntervalAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2IntervalPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2Interval ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_INTERVAL);
   } else {
      atp = AsnLinkType(orig, ID2_INTERVAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2IntervalNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_INTERVAL_gi) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gi = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_INTERVAL_start) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> start = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_INTERVAL_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> length = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2IntervalFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2IntervalAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2IntervalAsnWrite(ID2IntervalPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_INTERVAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> gi;
   retval = AsnWrite(aip, ID2_INTERVAL_gi,  &av);
   av.intvalue = ptr -> start;
   retval = AsnWrite(aip, ID2_INTERVAL_start,  &av);
   av.intvalue = ptr -> length;
   retval = AsnWrite(aip, ID2_INTERVAL_length,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2PackedSeqIntsNew()
*
**************************************************/
NLM_EXTERN 
ID2PackedSeqIntsPtr LIBCALL
ID2PackedSeqIntsNew(void)
{
   ID2PackedSeqIntsPtr ptr = MemNew((size_t) sizeof(ID2PackedSeqInts));

   return ptr;

}


/**************************************************
*
*    ID2PackedSeqIntsFree()
*
**************************************************/
NLM_EXTERN 
ID2PackedSeqIntsPtr LIBCALL
ID2PackedSeqIntsFree(ID2PackedSeqIntsPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr -> ints, (AsnOptFreeFunc) ID2SeqRangeFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2PackedSeqIntsAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2PackedSeqIntsPtr LIBCALL
ID2PackedSeqIntsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2PackedSeqIntsPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2PackedSeqInts ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_PACKED_SEQ_INTS);
   } else {
      atp = AsnLinkType(orig, ID2_PACKED_SEQ_INTS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2PackedSeqIntsNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_PACKED_SEQ_INTS_gi) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gi = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_PACKED_SEQ_INTS_ints) {
      ptr -> ints = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ID2SeqRangeAsnRead, (AsnOptFreeFunc) ID2SeqRangeFree);
      if (isError && ptr -> ints == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2PackedSeqIntsFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2PackedSeqIntsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2PackedSeqIntsAsnWrite(ID2PackedSeqIntsPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_PACKED_SEQ_INTS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> gi;
   retval = AsnWrite(aip, ID2_PACKED_SEQ_INTS_gi,  &av);
   AsnGenericUserSeqOfAsnWrite(ptr -> ints, (AsnWriteFunc) ID2SeqRangeAsnWrite, aip, ID2_PACKED_SEQ_INTS_ints, ID2_PACKED_SEQ_INTS_ints_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2SeqRangeNew()
*
**************************************************/
NLM_EXTERN 
ID2SeqRangePtr LIBCALL
ID2SeqRangeNew(void)
{
   ID2SeqRangePtr ptr = MemNew((size_t) sizeof(ID2SeqRange));

   ptr -> length = 1;
   return ptr;

}


/**************************************************
*
*    ID2SeqRangeFree()
*
**************************************************/
NLM_EXTERN 
ID2SeqRangePtr LIBCALL
ID2SeqRangeFree(ID2SeqRangePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    ID2SeqRangeAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2SeqRangePtr LIBCALL
ID2SeqRangeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2SeqRangePtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2SeqRange ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_SEQ_RANGE);
   } else {
      atp = AsnLinkType(orig, ID2_SEQ_RANGE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2SeqRangeNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_SEQ_RANGE_start) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> start = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_SEQ_RANGE_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> length = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2SeqRangeFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2SeqRangeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2SeqRangeAsnWrite(ID2SeqRangePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_SEQ_RANGE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> start;
   retval = AsnWrite(aip, ID2_SEQ_RANGE_start,  &av);
   av.intvalue = ptr -> length;
   retval = AsnWrite(aip, ID2_SEQ_RANGE_length,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ID2ParamNew()
*
**************************************************/
NLM_EXTERN 
ID2ParamPtr LIBCALL
ID2ParamNew(void)
{
   ID2ParamPtr ptr = MemNew((size_t) sizeof(ID2Param));

   ptr -> type = 1;
   return ptr;

}


/**************************************************
*
*    ID2ParamFree()
*
**************************************************/
NLM_EXTERN 
ID2ParamPtr LIBCALL
ID2ParamFree(ID2ParamPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   AsnGenericBaseSeqOfFree(ptr -> value ,ASNCODE_PTRVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    ID2ParamAsnRead()
*
**************************************************/
NLM_EXTERN 
ID2ParamPtr LIBCALL
ID2ParamAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ID2ParamPtr ptr;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ID2Param ::= (self contained) */
      atp = AsnReadId(aip, amp, ID2_PARAM);
   } else {
      atp = AsnLinkType(orig, ID2_PARAM);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ID2ParamNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ID2_PARAM_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_PARAM_value) {
      ptr -> value = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> value == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ID2_PARAM_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ID2ParamFree(ptr);
   goto ret;
}



/**************************************************
*
*    ID2ParamAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ID2ParamAsnWrite(ID2ParamPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! id2genAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ID2_PARAM);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, ID2_PARAM_name,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> value ,ASNCODE_PTRVAL_SLOT, aip, ID2_PARAM_value, ID2_PARAM_value_E);
   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, ID2_PARAM_type,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

