#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objseq.h>
#include <objcdd.h>

static Boolean loaded = FALSE;

#include <cdd.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objcddAsnLoad(void)
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
*    Generated object loaders for Module NCBI-Cdd
*    Generated using ASNCODE Revision: 6.8 at Oct 19, 1999  1:13 PM
*
**************************************************/


/**************************************************
*
*    CddIdFree()
*
**************************************************/
NLM_EXTERN 
CddIdPtr LIBCALL
CddIdFree(ValNodePtr anp)
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
   case CddId_gid:
      GlobalIdFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    CddIdAsnRead()
*
**************************************************/
NLM_EXTERN 
CddIdPtr LIBCALL
CddIdAsnRead(AsnIoPtr aip, AsnTypePtr orig)
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
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CddId ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD_ID);
   } else {
      atp = AsnLinkType(orig, CDD_ID);    /* link in local tree */
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
   if (atp == CDD_ID_uid) {
      choice = CddId_uid;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == CDD_ID_gid) {
      choice = CddId_gid;
      func = (AsnReadFunc) GlobalIdAsnRead;
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
*    CddIdAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddIdAsnWrite(CddIdPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, CDD_ID);   /* link local tree */
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
   case CddId_uid:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, CDD_ID_uid, &av);
      break;
   case CddId_gid:
      writetype = CDD_ID_gid;
      func = (AsnWriteFunc) GlobalIdAsnWrite;
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
*    CddIdSetFree()
*
**************************************************/
NLM_EXTERN 
CddIdSetPtr LIBCALL
CddIdSetFree(CddIdSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) CddIdFree);
   return NULL;
}


/**************************************************
*
*    CddIdSetAsnRead()
*
**************************************************/
NLM_EXTERN 
CddIdSetPtr LIBCALL
CddIdSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CddIdSetPtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CddIdSet ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD_ID_SET);
   } else {
      atp = AsnLinkType(orig, CDD_ID_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) CddIdAsnRead, (AsnOptFreeFunc) CddIdFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CddIdSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    CddIdSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddIdSetAsnWrite(CddIdSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDD_ID_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) CddIdAsnWrite, aip, atp, CDD_ID_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CddNew()
*
**************************************************/
NLM_EXTERN 
CddPtr LIBCALL
CddNew(void)
{
   CddPtr ptr = MemNew((size_t) sizeof(Cdd));

   return ptr;

}


/**************************************************
*
*    CddFree()
*
**************************************************/
NLM_EXTERN 
CddPtr LIBCALL
CddFree(CddPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   AsnGenericBaseSeqOfFree(ptr -> othernames ,ASNCODE_PTRVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> categories ,ASNCODE_PTRVAL_SLOT);
   CddIdFree(ptr -> id);
   AsnGenericChoiceSeqOfFree(ptr -> description, (AsnOptFreeFunc) CddDescrFree);
   DateFree(ptr -> create_date);
   MemFree(ptr -> source);
   SeqAnnotFree(ptr -> seqannot);
   BiostrucFeatureSetFree(ptr -> strucalign);
   return MemFree(ptr);
}


/**************************************************
*
*    CddAsnRead()
*
**************************************************/
NLM_EXTERN 
CddPtr LIBCALL
CddAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CddPtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Cdd ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD);
   } else {
      atp = AsnLinkType(orig, CDD);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CddNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CDD_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_othernames) {
      ptr -> othernames = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> othernames == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_categories) {
      ptr -> categories = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> categories == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_id) {
      ptr -> id = CddIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_description) {
      ptr -> description = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) CddDescrAsnRead, (AsnOptFreeFunc) CddDescrFree);
      if (isError && ptr -> description == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_create_date) {
      ptr -> create_date = DateAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_source) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> source = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_seqannot) {
      ptr -> seqannot = SeqAnnotAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_strucalign) {
      ptr -> strucalign = BiostrucFeatureSetAsnRead(aip, atp);
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
   ptr = CddFree(ptr);
   goto ret;
}



/**************************************************
*
*    CddAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddAsnWrite(CddPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDD);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, CDD_name,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> othernames ,ASNCODE_PTRVAL_SLOT, aip, CDD_othernames, CDD_othernames_E);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> categories ,ASNCODE_PTRVAL_SLOT, aip, CDD_categories, CDD_categories_E);
   if (ptr -> id != NULL) {
      if ( ! CddIdAsnWrite(ptr -> id, aip, CDD_id)) {
         goto erret;
      }
   }
   AsnGenericChoiceSeqOfAsnWrite(ptr -> description, (AsnWriteFunc) CddDescrAsnWrite, aip, CDD_description, CDD_description_E);
   if (ptr -> create_date != NULL) {
      if ( ! DateAsnWrite(ptr -> create_date, aip, CDD_create_date)) {
         goto erret;
      }
   }
   if (ptr -> source != NULL) {
      av.ptrvalue = ptr -> source;
      retval = AsnWrite(aip, CDD_source,  &av);
   }
   if (ptr -> seqannot != NULL) {
      if ( ! SeqAnnotAsnWrite(ptr -> seqannot, aip, CDD_seqannot)) {
         goto erret;
      }
   }
   if (ptr -> strucalign != NULL) {
      if ( ! BiostrucFeatureSetAsnWrite(ptr -> strucalign, aip, CDD_strucalign)) {
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
*    CddSetFree()
*
**************************************************/
NLM_EXTERN 
CddSetPtr LIBCALL
CddSetFree(CddSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) CddFree);
   return NULL;
}


/**************************************************
*
*    CddSetAsnRead()
*
**************************************************/
NLM_EXTERN 
CddSetPtr LIBCALL
CddSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CddSetPtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CddSet ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD_SET);
   } else {
      atp = AsnLinkType(orig, CDD_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) CddAsnRead, (AsnOptFreeFunc) CddFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CddSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    CddSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddSetAsnWrite(CddSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDD_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) CddAsnWrite, aip, atp, CDD_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CddTreeNew()
*
**************************************************/
NLM_EXTERN 
CddTreePtr LIBCALL
CddTreeNew(void)
{
   CddTreePtr ptr = MemNew((size_t) sizeof(CddTree));

   return ptr;

}


/**************************************************
*
*    CddTreeFree()
*
**************************************************/
NLM_EXTERN 
CddTreePtr LIBCALL
CddTreeFree(CddTreePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   CddIdFree(ptr -> id);
   CddIdSetFree(ptr -> parents);
   CddIdSetFree(ptr -> children);
   MemFree(ptr -> name);
   return MemFree(ptr);
}


/**************************************************
*
*    CddTreeAsnRead()
*
**************************************************/
NLM_EXTERN 
CddTreePtr LIBCALL
CddTreeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CddTreePtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CddTree ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD_TREE);
   } else {
      atp = AsnLinkType(orig, CDD_TREE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CddTreeNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CDD_TREE_id) {
      ptr -> id = CddIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_TREE_parents) {
      ptr -> parents = CddIdSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_TREE_children) {
      ptr -> children = CddIdSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_TREE_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
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
   ptr = CddTreeFree(ptr);
   goto ret;
}



/**************************************************
*
*    CddTreeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddTreeAsnWrite(CddTreePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDD_TREE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> id != NULL) {
      if ( ! CddIdAsnWrite(ptr -> id, aip, CDD_TREE_id)) {
         goto erret;
      }
   }
   if (ptr -> parents != NULL) {
      if ( ! CddIdSetAsnWrite(ptr -> parents, aip, CDD_TREE_parents)) {
         goto erret;
      }
   }
   if (ptr -> children != NULL) {
      if ( ! CddIdSetAsnWrite(ptr -> children, aip, CDD_TREE_children)) {
         goto erret;
      }
   }
   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, CDD_TREE_name,  &av);
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
*    CddTreeSetFree()
*
**************************************************/
NLM_EXTERN 
CddTreeSetPtr LIBCALL
CddTreeSetFree(CddTreeSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) CddTreeFree);
   return NULL;
}


/**************************************************
*
*    CddTreeSetAsnRead()
*
**************************************************/
NLM_EXTERN 
CddTreeSetPtr LIBCALL
CddTreeSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CddTreeSetPtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CddTreeSet ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD_TREE_SET);
   } else {
      atp = AsnLinkType(orig, CDD_TREE_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) CddTreeAsnRead, (AsnOptFreeFunc) CddTreeFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CddTreeSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    CddTreeSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddTreeSetAsnWrite(CddTreeSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CDD_TREE_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) CddTreeAsnWrite, aip, atp, CDD_TREE_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GlobalIdNew()
*
**************************************************/
NLM_EXTERN 
GlobalIdPtr LIBCALL
GlobalIdNew(void)
{
   GlobalIdPtr ptr = MemNew((size_t) sizeof(GlobalId));

   return ptr;

}


/**************************************************
*
*    GlobalIdFree()
*
**************************************************/
NLM_EXTERN 
GlobalIdPtr LIBCALL
GlobalIdFree(GlobalIdPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> accession);
   MemFree(ptr -> release);
   MemFree(ptr -> database);
   return MemFree(ptr);
}


/**************************************************
*
*    GlobalIdAsnRead()
*
**************************************************/
NLM_EXTERN 
GlobalIdPtr LIBCALL
GlobalIdAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GlobalIdPtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GlobalId ::= (self contained) */
      atp = AsnReadId(aip, amp, GLOBAL_ID);
   } else {
      atp = AsnLinkType(orig, GLOBAL_ID);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GlobalIdNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GLOBAL_ID_accession) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> accession = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GLOBAL_ID_release) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> release = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GLOBAL_ID_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> version = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GLOBAL_ID_database) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> database = av.ptrvalue;
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
   ptr = GlobalIdFree(ptr);
   goto ret;
}



/**************************************************
*
*    GlobalIdAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GlobalIdAsnWrite(GlobalIdPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GLOBAL_ID);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> accession != NULL) {
      av.ptrvalue = ptr -> accession;
      retval = AsnWrite(aip, GLOBAL_ID_accession,  &av);
   }
   if (ptr -> release != NULL) {
      av.ptrvalue = ptr -> release;
      retval = AsnWrite(aip, GLOBAL_ID_release,  &av);
   }
   av.intvalue = ptr -> version;
   retval = AsnWrite(aip, GLOBAL_ID_version,  &av);
   if (ptr -> database != NULL) {
      av.ptrvalue = ptr -> database;
      retval = AsnWrite(aip, GLOBAL_ID_database,  &av);
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
*    CddDescrFree()
*
**************************************************/
NLM_EXTERN 
CddDescrPtr LIBCALL
CddDescrFree(ValNodePtr anp)
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
   case CddDescr_comment:
      MemFree(anp -> data.ptrvalue);
      break;
   case CddDescr_reference:
      PubFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    CddDescrAsnRead()
*
**************************************************/
NLM_EXTERN 
CddDescrPtr LIBCALL
CddDescrAsnRead(AsnIoPtr aip, AsnTypePtr orig)
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
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CddDescr ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD_DESCR);
   } else {
      atp = AsnLinkType(orig, CDD_DESCR);    /* link in local tree */
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
   if (atp == CDD_DESCR_comment) {
      choice = CddDescr_comment;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == CDD_DESCR_reference) {
      choice = CddDescr_reference;
      func = (AsnReadFunc) PubAsnRead;
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
*    CddDescrAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddDescrAsnWrite(CddDescrPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objcddAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, CDD_DESCR);   /* link local tree */
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
   case CddDescr_comment:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, CDD_DESCR_comment, &av);
      break;
   case CddDescr_reference:
      writetype = CDD_DESCR_reference;
      func = (AsnWriteFunc) PubAsnWrite;
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
