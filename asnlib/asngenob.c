/*  asngenob.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* File Name:  asngenob.c
*
* Author:  Sirotkin/Epstein
*
* Version Creation Date: 4/21/94
*
* $Revision: 6.1 $
*
* File Description:
*   Generic routines shared by object loaders which are automatically
*   generated by ASNCODE.
*   
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: asngenob.c,v $
* Revision 6.1  1997/10/29 02:40:56  vakatov
* Type castings to pass through the C++ compiler
*
* Revision 6.0  1997/08/25 18:09:57  madden
* Revision changed to 6.0
*
* Revision 5.1  1996/12/03 21:43:48  vakatov
* Adopted for 32-bit MS-Windows DLLs
*
 * Revision 5.0  1996/05/28  14:00:29  ostell
 * Set to revision 5.0
 *
 * Revision 4.0  1995/07/26  13:47:38  ostell
 * force revision to 4.0
 *
 * Revision 1.3  1995/05/15  18:38:28  ostell
 * added Log line
 *
*
* ==========================================================================
*/

#include <asn.h>

typedef struct struct__Generic_linked_list {
    struct struct__Generic_linked_list PNTR next;
} AsnGenericLinkedList, PNTR AsnGenericLinkListPtr;


NLM_EXTERN ValNodePtr LIBCALL AsnGenericBaseSeqOfAsnRead (AsnIoPtr aip, AsnModulePtr amp, AsnTypePtr orig, int whichvalslot, BoolPtr isError)
{
   AsnTypePtr      start_atp;
   ValNodePtr      current;
   ValNodePtr      head = NULL;
   ValNodePtr      prev = NULL;
   AsnTypePtr      atp = orig;
   DataVal         av;

   if (isError != NULL)
      *isError = FALSE;

   if (aip == NULL)
      return NULL;

   if (AsnReadVal (aip, atp, &av) <= 0)	/* START_STRUCT */
      goto erret;

   start_atp = atp;

   while ((atp = AsnReadId (aip, amp, atp)) != start_atp) {
      if (atp == NULL)
	 goto erret;

      current = ValNodeNew (prev);
      if (current == NULL)
	 goto erret;
      switch (whichvalslot) {
      case ASNCODE_PTRVAL_SLOT:
      case ASNCODE_BYTEVAL_SLOT:
	 if (AsnReadVal (aip, atp, &av) <= 0)
	    goto erret;
	 current->data.ptrvalue = av.ptrvalue;
	 break;
      case ASNCODE_REALVAL_SLOT:
	 if (AsnReadVal (aip, atp, &av) <= 0)
	    goto erret;
	 current->data.realvalue = av.realvalue;
	 break;
      case ASNCODE_INTVAL_SLOT:
	 if (AsnReadVal (aip, atp, &av) <= 0)
	    goto erret;
	 current->data.intvalue = av.intvalue;
	 break;
      case ASNCODE_BOOLVAL_SLOT:
	 if (AsnReadVal (aip, atp, &av) <= 0)
	    goto erret;
	 current->data.boolvalue = av.boolvalue;
	 break;
      }

      if (head == NULL)
	 head = current;
      else
	 prev->next = current;
      prev = current;
   }
   if (AsnReadVal (aip, atp, &av) <= 0)	/* END_STRUCT */
      goto erret;

ret:
   return head;

erret:
   head = ValNodeFree (head);
   if (isError != NULL)
      *isError = TRUE;
   goto ret;
}

NLM_EXTERN Boolean LIBCALL AsnGenericBaseSeqOfAsnWrite (ValNodePtr ptr, int whichvalslot, AsnIoPtr aip, AsnTypePtr bag_atp, AsnTypePtr element_atp)
{
   Boolean         retval = FALSE;

   if (ptr == NULL && bag_atp->optional)
      return TRUE;

   if (!AsnOpenStruct (aip, bag_atp, ptr))
      goto ret;

   if (ptr != NULL) {
      ValNodePtr      current;

      for (current = ptr; current; current = current->next) {
	 switch (whichvalslot) {
	 case ASNCODE_PTRVAL_SLOT:
	 case ASNCODE_BYTEVAL_SLOT:
	 case ASNCODE_REALVAL_SLOT:
	 case ASNCODE_INTVAL_SLOT:
	 case ASNCODE_BOOLVAL_SLOT:
	    if (AsnWrite (aip, element_atp, &(current->data)) <= 0)
	       goto ret;
	    break;
	 }
      }

   }
   if (!AsnCloseStruct (aip, bag_atp, ptr)) {
      goto ret;
   }
   retval = TRUE;

ret:
   return retval;
}


NLM_EXTERN Boolean LIBCALL AsnGenericBaseSeqOfFree (ValNodePtr ptr, int whichvalslot)
{
   Boolean         retval = TRUE;
   ValNodePtr      current, nextcur;

   for (current = ptr; current; current = nextcur) {
      nextcur = current->next;
      switch (whichvalslot) {
      case ASNCODE_BYTEVAL_SLOT:
	 BSFree ((ByteStorePtr)current->data.ptrvalue);
	 break;
      case ASNCODE_PTRVAL_SLOT:
	 MemFree (current->data.ptrvalue);
	 current->data.ptrvalue = NULL;
	 break;
      case ASNCODE_REALVAL_SLOT:
      case ASNCODE_INTVAL_SLOT:
      case ASNCODE_BOOLVAL_SLOT:
	 /* No-op */
	 break;
      }
   }

   ValNodeFree (ptr);
   return retval;
}


NLM_EXTERN Pointer AsnGenericUserSeqOfAsnRead (AsnIoPtr aip, AsnModulePtr amp, AsnTypePtr orig, BoolPtr isError, AsnReadFunc readfunc, AsnOptFreeFunc freefunc)
{
   AsnTypePtr      start_atp;
   AsnGenericLinkListPtr current;
   AsnGenericLinkListPtr head = NULL;
   AsnGenericLinkListPtr prev = NULL;
   AsnTypePtr      atp = orig;
   DataVal         av;

   if (isError != NULL)
      *isError = FALSE;

   if (aip == NULL)
      return NULL;

   if (AsnReadVal (aip, atp, &av) <= 0)	/* START_STRUCT */
      goto erret;

   start_atp = atp;

   while ((atp = AsnReadId (aip, amp, atp)) != start_atp) {
      if (atp == NULL)
	 goto erret;

      current = (AsnGenericLinkListPtr) readfunc (aip, atp);
      if (current == NULL)
	 goto erret;

      if (head == NULL)
	 head = current;
      else
	 prev->next = current;
      prev = current;
   }
   if (AsnReadVal (aip, atp, &av) <= 0)	/* END_STRUCT */
      goto erret;

ret:
   return (Pointer) head;

erret:
   head = (AsnGenericLinkListPtr) freefunc (head);
   if (isError != NULL)
      *isError = TRUE;

   goto ret;
}


NLM_EXTERN Boolean LIBCALL AsnGenericUserSeqOfAsnWrite (Pointer ptr, AsnWriteFunc writefunc, AsnIoPtr aip, AsnTypePtr bag_atp, AsnTypePtr element_atp)
{
   Boolean         retval = FALSE;

   if (ptr == NULL && bag_atp->optional)
      return TRUE;

   if (!AsnOpenStruct (aip, bag_atp, ptr))
      goto ret;

   if (ptr != NULL) {
      AsnGenericLinkListPtr current;

      for (current = (AsnGenericLinkListPtr) ptr; current; current = current->next) {
	 if (!writefunc (current, aip, element_atp))
	    goto ret;
      }

   }
   if (!AsnCloseStruct (aip, bag_atp, ptr)) {
      goto ret;
   }
   retval = TRUE;

ret:
   return retval;
}


NLM_EXTERN Boolean LIBCALL AsnGenericUserSeqOfFree (Pointer ptr, AsnOptFreeFunc freefunc)
{
   Boolean         retval = TRUE;
   AsnGenericLinkListPtr current, nextcur;

   for (current = (AsnGenericLinkListPtr) ptr; current; current = nextcur) {
      nextcur = current->next;
      freefunc (current);
   }

   return retval;
}


NLM_EXTERN Pointer LIBCALL AsnGenericChoiceSeqOfAsnRead (AsnIoPtr aip, AsnModulePtr amp, AsnTypePtr orig, BoolPtr isError, AsnReadFunc readfunc, AsnOptFreeFunc freefunc)
{
   AsnTypePtr      start_atp;
   ValNodePtr      current;
   ValNodePtr      head = NULL;
   ValNodePtr      prev = NULL;
   AsnTypePtr      atp = orig;
   DataVal         av;

   if (isError != NULL)
      *isError = FALSE;

   if (aip == NULL)
      return NULL;

   if (AsnReadVal (aip, atp, &av) <= 0)	/* START_STRUCT */
      goto erret;

   start_atp = atp;

   while ((atp = AsnReadId (aip, amp, atp)) != start_atp) {
      if (atp == NULL)
	 goto erret;

      current = (ValNodePtr) readfunc (aip, atp);
      if (current == NULL)
	 goto erret;

      if (head == NULL)
	 head = current;
      else
	 prev->next = current;
      prev = current;
   }
   if (AsnReadVal (aip, atp, &av) <= 0)	/* END_STRUCT */
      goto erret;

ret:
   return (Pointer) head;

erret:
   head = (ValNodePtr) freefunc (head);
   if (isError != NULL)
      *isError = TRUE;

   goto ret;
}


NLM_EXTERN Boolean LIBCALL AsnGenericChoiceSeqOfAsnWrite (Pointer ptr, AsnWriteFunc writefunc, AsnIoPtr aip, AsnTypePtr bag_atp, AsnTypePtr element_atp)
{
   Boolean         retval = FALSE;

   if (ptr == NULL && bag_atp->optional)
      return TRUE;

   if (!AsnOpenStruct (aip, bag_atp, ptr))
      goto ret;

   if (ptr != NULL) {
      ValNodePtr      current;

      for (current = (ValNodePtr) ptr; current; current = current->next) {
	 if (!writefunc (current, aip, element_atp))
	    goto ret;
      }

   }
   if (!AsnCloseStruct (aip, bag_atp, ptr)) {
      goto ret;
   }
   retval = TRUE;

ret:
   return retval;
}

NLM_EXTERN Boolean LIBCALL AsnGenericChoiceSeqOfFree (Pointer ptr, AsnOptFreeFunc freefunc)
{
   Boolean         retval = TRUE;
   ValNodePtr      current, nextcur;

   for (current = (ValNodePtr) ptr; current; current = nextcur) {
      nextcur = current->next;
      freefunc (current);
   }

   return retval;
}
