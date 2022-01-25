#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objsset.h>
#include <objmmdb1.h>
#include <objmmdb3.h>
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
*    Generated using ASNCODE Revision: 6.9 at May 12, 2000  4:52 PM
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
   CddIdSetFree(ptr -> id);
   CddDescrSetFree(ptr -> description);
   AsnGenericUserSeqOfFree(ptr -> seqannot, (AsnOptFreeFunc) SeqAnnotFree);
   BiostrucAnnotSetFree(ptr -> features);
   SeqEntryFree(ptr -> sequences);
   SeqIntFree(ptr -> profile_range);
   BioseqFree(ptr -> trunc_master);
   MatrixFree(ptr -> posfreq);
   MatrixFree(ptr -> scoremat);
   TriangleFree(ptr -> distance);
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
   if (atp == CDD_id) {
      ptr -> id = CddIdSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_description) {
      ptr -> description = CddDescrSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_seqannot) {
      ptr -> seqannot = SeqAnnotSetAsnRead(aip,CDD_seqannot,CDD_seqannot_E);
      if (isError && ptr -> seqannot == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_features) {
      ptr -> features = BiostrucAnnotSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_sequences) {
      ptr -> sequences = SeqEntryAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_profile_range) {
      ptr -> profile_range = SeqIntAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_trunc_master) {
      ptr -> trunc_master = BioseqAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_posfreq) {
      ptr -> posfreq = MatrixAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_scoremat) {
      ptr -> scoremat = MatrixAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_distance) {
      ptr -> distance = TriangleAsnRead(aip, atp);
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
   if (ptr -> id != NULL) {
      if ( ! CddIdSetAsnWrite(ptr -> id, aip, CDD_id)) {
         goto erret;
      }
   }
   if (ptr -> description != NULL) {
      if ( ! CddDescrSetAsnWrite(ptr -> description, aip, CDD_description)) {
         goto erret;
      }
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> seqannot, (AsnWriteFunc) SeqAnnotAsnWrite, aip, CDD_seqannot, CDD_seqannot_E);
   if (ptr -> features != NULL) {
      if ( ! BiostrucAnnotSetAsnWrite(ptr -> features, aip, CDD_features)) {
         goto erret;
      }
   }
   if (ptr -> sequences != NULL) {
      if ( ! SeqEntryAsnWrite(ptr -> sequences, aip, CDD_sequences)) {
         goto erret;
      }
   }
   if (ptr -> profile_range != NULL) {
      if ( ! SeqIntAsnWrite(ptr -> profile_range, aip, CDD_profile_range)) {
         goto erret;
      }
   }
   if (ptr -> trunc_master != NULL) {
      if ( ! BioseqAsnWrite(ptr -> trunc_master, aip, CDD_trunc_master)) {
         goto erret;
      }
   }
   if (ptr -> posfreq != NULL) {
      if ( ! MatrixAsnWrite(ptr -> posfreq, aip, CDD_posfreq)) {
         goto erret;
      }
   }
   if (ptr -> scoremat != NULL) {
      if ( ! MatrixAsnWrite(ptr -> scoremat, aip, CDD_scoremat)) {
         goto erret;
      }
   }
   if (ptr -> distance != NULL) {
      if ( ! TriangleAsnWrite(ptr -> distance, aip, CDD_distance)) {
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
   MemFree(ptr -> name);
   CddIdSetFree(ptr -> id);
   CddDescrSetFree(ptr -> description);
   CddIdSetFree(ptr -> parents);
   CddIdSetFree(ptr -> children);
   CddIdSetFree(ptr -> siblings);
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

   if (atp == CDD_TREE_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_TREE_id) {
      ptr -> id = CddIdSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CDD_TREE_description) {
      ptr -> description = CddDescrSetAsnRead(aip, atp);
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
   if (atp == CDD_TREE_siblings) {
      ptr -> siblings = CddIdSetAsnRead(aip, atp);
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

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, CDD_TREE_name,  &av);
   }
   if (ptr -> id != NULL) {
      if ( ! CddIdSetAsnWrite(ptr -> id, aip, CDD_TREE_id)) {
         goto erret;
      }
   }
   if (ptr -> description != NULL) {
      if ( ! CddDescrSetAsnWrite(ptr -> description, aip, CDD_TREE_description)) {
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
   if (ptr -> siblings != NULL) {
      if ( ! CddIdSetAsnWrite(ptr -> siblings, aip, CDD_TREE_siblings)) {
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
   case CddDescr_othername:
      MemFree(anp -> data.ptrvalue);
      break;
   case CddDescr_category:
      MemFree(anp -> data.ptrvalue);
      break;
   case CddDescr_comment:
      MemFree(anp -> data.ptrvalue);
      break;
   case CddDescr_reference:
      PubFree(anp -> data.ptrvalue);
      break;
   case CddDescr_create_date:
      DateFree(anp -> data.ptrvalue);
      break;
   case CddDescr_tax_source:
      OrgRefFree(anp -> data.ptrvalue);
      break;
   case CddDescr_source:
      MemFree(anp -> data.ptrvalue);
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
   if (atp == CDD_DESCR_othername) {
      choice = CddDescr_othername;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == CDD_DESCR_category) {
      choice = CddDescr_category;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == CDD_DESCR_comment) {
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
   else if (atp == CDD_DESCR_create_date) {
      choice = CddDescr_create_date;
      func = (AsnReadFunc) DateAsnRead;
   }
   else if (atp == CDD_DESCR_tax_source) {
      choice = CddDescr_tax_source;
      func = (AsnReadFunc) OrgRefAsnRead;
   }
   else if (atp == CDD_DESCR_source) {
      choice = CddDescr_source;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == CDD_DESCR_status) {
      choice = CddDescr_status;
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
   case CddDescr_othername:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, CDD_DESCR_othername, &av);
      break;
   case CddDescr_category:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, CDD_DESCR_category, &av);
      break;
   case CddDescr_comment:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, CDD_DESCR_comment, &av);
      break;
   case CddDescr_reference:
      writetype = CDD_DESCR_reference;
      func = (AsnWriteFunc) PubAsnWrite;
      break;
   case CddDescr_create_date:
      writetype = CDD_DESCR_create_date;
      func = (AsnWriteFunc) DateAsnWrite;
      break;
   case CddDescr_tax_source:
      writetype = CDD_DESCR_tax_source;
      func = (AsnWriteFunc) OrgRefAsnWrite;
      break;
   case CddDescr_source:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, CDD_DESCR_source, &av);
      break;
   case CddDescr_status:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, CDD_DESCR_status, &av);
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
*    CddDescrSetFree()
*
**************************************************/
NLM_EXTERN 
CddDescrSetPtr LIBCALL
CddDescrSetFree(CddDescrSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr, (AsnOptFreeFunc) CddDescrFree);
   return NULL;
}


/**************************************************
*
*    CddDescrSetAsnRead()
*
**************************************************/
NLM_EXTERN 
CddDescrSetPtr LIBCALL
CddDescrSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CddDescrSetPtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CddDescrSet ::= (self contained) */
      atp = AsnReadId(aip, amp, CDD_DESCR_SET);
   } else {
      atp = AsnLinkType(orig, CDD_DESCR_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) CddDescrAsnRead, (AsnOptFreeFunc) CddDescrFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CddDescrSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    CddDescrSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CddDescrSetAsnWrite(CddDescrSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
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

   atp = AsnLinkType(orig, CDD_DESCR_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericChoiceSeqOfAsnWrite(ptr , (AsnWriteFunc) CddDescrAsnWrite, aip, atp, CDD_DESCR_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    MatrixNew()
*
**************************************************/
NLM_EXTERN 
MatrixPtr LIBCALL
MatrixNew(void)
{
   MatrixPtr ptr = MemNew((size_t) sizeof(Matrix));

   return ptr;

}


/**************************************************
*
*    MatrixFree()
*
**************************************************/
NLM_EXTERN 
MatrixPtr LIBCALL
MatrixFree(MatrixPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> row_labels ,ASNCODE_PTRVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> columns ,ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    MatrixAsnRead()
*
**************************************************/
NLM_EXTERN 
MatrixPtr LIBCALL
MatrixAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   MatrixPtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Matrix ::= (self contained) */
      atp = AsnReadId(aip, amp, MATRIX);
   } else {
      atp = AsnLinkType(orig, MATRIX);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = MatrixNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MATRIX_ncolumns) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> ncolumns = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MATRIX_nrows) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> nrows = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MATRIX_row_labels) {
      ptr -> row_labels = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> row_labels == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MATRIX_scale_factor) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> scale_factor = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MATRIX_columns) {
      ptr -> columns = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> columns == NULL) {
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
   ptr = MatrixFree(ptr);
   goto ret;
}



/**************************************************
*
*    MatrixAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
MatrixAsnWrite(MatrixPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
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

   atp = AsnLinkType(orig, MATRIX);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> ncolumns;
   retval = AsnWrite(aip, MATRIX_ncolumns,  &av);
   av.intvalue = ptr -> nrows;
   retval = AsnWrite(aip, MATRIX_nrows,  &av);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> row_labels ,ASNCODE_PTRVAL_SLOT, aip, MATRIX_row_labels, MATRIX_row_labels_E);
   av.intvalue = ptr -> scale_factor;
   retval = AsnWrite(aip, MATRIX_scale_factor,  &av);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> columns ,ASNCODE_INTVAL_SLOT, aip, MATRIX_columns, MATRIX_columns_E);
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
*    TriangleNew()
*
**************************************************/
NLM_EXTERN 
TrianglePtr LIBCALL
TriangleNew(void)
{
   TrianglePtr ptr = MemNew((size_t) sizeof(Triangle));

   return ptr;

}


/**************************************************
*
*    TriangleFree()
*
**************************************************/
NLM_EXTERN 
TrianglePtr LIBCALL
TriangleFree(TrianglePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ScoreSetFree(ptr -> scores);
   return MemFree(ptr);
}


/**************************************************
*
*    TriangleAsnRead()
*
**************************************************/
NLM_EXTERN 
TrianglePtr LIBCALL
TriangleAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TrianglePtr ptr;

   if (! loaded)
   {
      if (! objcddAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Triangle ::= (self contained) */
      atp = AsnReadId(aip, amp, TRIANGLE);
   } else {
      atp = AsnLinkType(orig, TRIANGLE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = TriangleNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == TRIANGLE_nelements) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> nelements = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TRIANGLE_scores) {
      ptr -> scores = ScoreSetAsnRead(aip, atp);
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
   ptr = TriangleFree(ptr);
   goto ret;
}



/**************************************************
*
*    TriangleAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TriangleAsnWrite(TrianglePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
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

   atp = AsnLinkType(orig, TRIANGLE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> nelements;
   retval = AsnWrite(aip, TRIANGLE_nelements,  &av);
   if (ptr -> scores != NULL) {
      if ( ! ScoreSetAsnWrite(ptr -> scores, aip, TRIANGLE_scores)) {
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

