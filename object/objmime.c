/*  objmime.c
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
* File Name:  objmime.c
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: objmime.c,v $
* Revision 6.4  1998/12/07 16:29:28  ywang
* add object loaded for mime type Biostruc-seqs
*
* ==========================================================================
*/

#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

/* hand change add include --- lyg */
#include <mapmime.h>
#include <objmime.h>

static Boolean loaded = FALSE;

/*hand change from ncbimime.h -- lyg & yanli */
#include <asnmime.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objmimeAsnLoad(void)
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
*    Generated object loaders for Module NCBI-Mime
*    Generated using ASNCODE Revision: 6.5 at Dec 4, 1998  2:11 PM
*
**************************************************/


/**************************************************
*
*    NcbiMimeAsn1Free()
*
**************************************************/
NLM_EXTERN 
NcbiMimeAsn1Ptr LIBCALL
NcbiMimeAsn1Free(ValNodePtr anp)
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
   case NcbiMimeAsn1_entrez:
      EntrezGeneralFree(anp -> data.ptrvalue);
      break;
   case NcbiMimeAsn1_alignstruc:
      BiostrucAlignFree(anp -> data.ptrvalue);
      break;
   case NcbiMimeAsn1_alignseq:
      BiostrucAlignSeqFree(anp -> data.ptrvalue);
      break;
   case NcbiMimeAsn1_strucseq:
      BiostrucSeqFree(anp -> data.ptrvalue);
      break;
   case NcbiMimeAsn1_strucseqs:
      BiostrucSeqsFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    NcbiMimeAsn1AsnRead()
*
**************************************************/
NLM_EXTERN 
NcbiMimeAsn1Ptr LIBCALL
NcbiMimeAsn1AsnRead(AsnIoPtr aip, AsnTypePtr orig)
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
      if (! objmimeAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* NcbiMimeAsn1 ::= (self contained) */
      atp = AsnReadId(aip, amp, NCBI_MIME_ASN1);
   } else {
      atp = AsnLinkType(orig, NCBI_MIME_ASN1);    /* link in local tree */
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
   if (atp == NCBI_MIME_ASN1_entrez) {
      choice = NcbiMimeAsn1_entrez;
      func = (AsnReadFunc) EntrezGeneralAsnRead;
   }
   else if (atp == NCBI_MIME_ASN1_alignstruc) {
      choice = NcbiMimeAsn1_alignstruc;
      func = (AsnReadFunc) BiostrucAlignAsnRead;
   }
   else if (atp == NCBI_MIME_ASN1_alignseq) {
      choice = NcbiMimeAsn1_alignseq;
      func = (AsnReadFunc) BiostrucAlignSeqAsnRead;
   }
   else if (atp == NCBI_MIME_ASN1_strucseq) {
      choice = NcbiMimeAsn1_strucseq;
      func = (AsnReadFunc) BiostrucSeqAsnRead;
   }
   else if (atp == NCBI_MIME_ASN1_strucseqs) {
      choice = NcbiMimeAsn1_strucseqs;
      func = (AsnReadFunc) BiostrucSeqsAsnRead;
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
*    NcbiMimeAsn1AsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
NcbiMimeAsn1AsnWrite(NcbiMimeAsn1Ptr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmimeAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, NCBI_MIME_ASN1);   /* link local tree */
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
   case NcbiMimeAsn1_entrez:
      writetype = NCBI_MIME_ASN1_entrez;
      func = (AsnWriteFunc) EntrezGeneralAsnWrite;
      break;
   case NcbiMimeAsn1_alignstruc:
      writetype = NCBI_MIME_ASN1_alignstruc;
      func = (AsnWriteFunc) BiostrucAlignAsnWrite;
      break;
   case NcbiMimeAsn1_alignseq:
      writetype = NCBI_MIME_ASN1_alignseq;
      func = (AsnWriteFunc) BiostrucAlignSeqAsnWrite;
      break;
   case NcbiMimeAsn1_strucseq:
      writetype = NCBI_MIME_ASN1_strucseq;
      func = (AsnWriteFunc) BiostrucSeqAsnWrite;
      break;
   case NcbiMimeAsn1_strucseqs:
      writetype = NCBI_MIME_ASN1_strucseqs;
      func = (AsnWriteFunc) BiostrucSeqsAsnWrite;
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
*    EntrezGeneralNew()
*
**************************************************/
NLM_EXTERN 
EntrezGeneralPtr LIBCALL
EntrezGeneralNew(void)
{
   EntrezGeneralPtr ptr = MemNew((size_t) sizeof(EntrezGeneral));

   return ptr;

}


/**************************************************
*
*    EntrezGeneralFree()
*
**************************************************/
NLM_EXTERN 
EntrezGeneralPtr LIBCALL
EntrezGeneralFree(EntrezGeneralPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> title);
   Data_dataFree(ptr -> Data_data);
   MemFree(ptr -> location);
   return MemFree(ptr);
}


/**************************************************
*
*    Data_dataFree()
*
**************************************************/
static 
Data_dataPtr LIBCALL
Data_dataFree(ValNodePtr anp)
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
   case Data_data_ml:
      MedlineEntryFree(anp -> data.ptrvalue);
      break;
   case Data_data_prot:
      SeqEntryFree(anp -> data.ptrvalue);
      break;
   case Data_data_nuc:
      SeqEntryFree(anp -> data.ptrvalue);
      break;
   case Data_data_genome:
      SeqEntryFree(anp -> data.ptrvalue);
      break;
   case Data_data_structure:
      BiostrucFree(anp -> data.ptrvalue);
      break;
   case Data_data_strucAnnot:
      BiostrucAnnotSetFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    EntrezGeneralAsnRead()
*
**************************************************/
NLM_EXTERN 
EntrezGeneralPtr LIBCALL
EntrezGeneralAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   EntrezGeneralPtr ptr;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* EntrezGeneral ::= (self contained) */
      atp = AsnReadId(aip, amp, ENTREZ_GENERAL);
   } else {
      atp = AsnLinkType(orig, ENTREZ_GENERAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = EntrezGeneralNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == ENTREZ_GENERAL_title) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> title = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ENTREZ_GENERAL_data) {
      ptr -> Data_data = Data_dataAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ENTREZ_GENERAL_style) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> style = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == ENTREZ_GENERAL_location) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> location = av.ptrvalue;
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
   ptr = EntrezGeneralFree(ptr);
   goto ret;
}



/**************************************************
*
*    Data_dataAsnRead()
*
**************************************************/
static 
Data_dataPtr LIBCALL
Data_dataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
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
      if (! objmimeAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Data_data ::= (self contained) */
      atp = AsnReadId(aip, amp, ENTREZ_GENERAL_data);
   } else {
      atp = AsnLinkType(orig, ENTREZ_GENERAL_data);    /* link in local tree */
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
   if (atp == ENTREZ_GENERAL_data_ml) {
      choice = Data_data_ml;
      func = (AsnReadFunc) MedlineEntryAsnRead;
   }
   else if (atp == ENTREZ_GENERAL_data_prot) {
      choice = Data_data_prot;
      func = (AsnReadFunc) SeqEntryAsnRead;
   }
   else if (atp == ENTREZ_GENERAL_data_nuc) {
      choice = Data_data_nuc;
      func = (AsnReadFunc) SeqEntryAsnRead;
   }
   else if (atp == ENTREZ_GENERAL_data_genome) {
      choice = Data_data_genome;
      func = (AsnReadFunc) SeqEntryAsnRead;
   }
   else if (atp == ENTREZ_GENERAL_data_structure) {
      choice = Data_data_structure;
      func = (AsnReadFunc) BiostrucAsnRead;
   }
   else if (atp == ENTREZ_GENERAL_data_strucAnnot) {
      choice = Data_data_strucAnnot;
      func = (AsnReadFunc) BiostrucAnnotSetAsnRead;
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
*    EntrezGeneralAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
EntrezGeneralAsnWrite(EntrezGeneralPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, ENTREZ_GENERAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> title != NULL) {
      av.ptrvalue = ptr -> title;
      retval = AsnWrite(aip, ENTREZ_GENERAL_title,  &av);
   }
   if (ptr -> Data_data != NULL) {
      if ( ! Data_dataAsnWrite(ptr -> Data_data, aip, ENTREZ_GENERAL_data)) {
         goto erret;
      }
   }
   av.intvalue = ptr -> style;
   retval = AsnWrite(aip, ENTREZ_GENERAL_style,  &av);
   if (ptr -> location != NULL) {
      av.ptrvalue = ptr -> location;
      retval = AsnWrite(aip, ENTREZ_GENERAL_location,  &av);
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
*    Data_dataAsnWrite()
*
**************************************************/
static Boolean LIBCALL 
Data_dataAsnWrite(Data_dataPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmimeAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, ENTREZ_GENERAL_data);   /* link local tree */
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
   case Data_data_ml:
      writetype = ENTREZ_GENERAL_data_ml;
      func = (AsnWriteFunc) MedlineEntryAsnWrite;
      break;
   case Data_data_prot:
      writetype = ENTREZ_GENERAL_data_prot;
      func = (AsnWriteFunc) SeqEntryAsnWrite;
      break;
   case Data_data_nuc:
      writetype = ENTREZ_GENERAL_data_nuc;
      func = (AsnWriteFunc) SeqEntryAsnWrite;
      break;
   case Data_data_genome:
      writetype = ENTREZ_GENERAL_data_genome;
      func = (AsnWriteFunc) SeqEntryAsnWrite;
      break;
   case Data_data_structure:
      writetype = ENTREZ_GENERAL_data_structure;
      func = (AsnWriteFunc) BiostrucAsnWrite;
      break;
   case Data_data_strucAnnot:
      writetype = ENTREZ_GENERAL_data_strucAnnot;
      func = (AsnWriteFunc) BiostrucAnnotSetAsnWrite;
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
*    BiostrucAlignNew()
*
**************************************************/
NLM_EXTERN 
BiostrucAlignPtr LIBCALL
BiostrucAlignNew(void)
{
   BiostrucAlignPtr ptr = MemNew((size_t) sizeof(BiostrucAlign));

   return ptr;

}


/**************************************************
*
*    BiostrucAlignFree()
*
**************************************************/
NLM_EXTERN 
BiostrucAlignPtr LIBCALL
BiostrucAlignFree(BiostrucAlignPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   BiostrucFree(ptr -> master);
   AsnGenericUserSeqOfFree(ptr -> slaves, (AsnOptFreeFunc) BiostrucFree);
   BiostrucAnnotSetFree(ptr -> alignments);
   AsnGenericChoiceSeqOfFree(ptr -> sequences, (AsnOptFreeFunc) SeqEntryFree);
   AsnGenericUserSeqOfFree(ptr -> seqalign, (AsnOptFreeFunc) SeqAnnotFree);
   return MemFree(ptr);
}


/**************************************************
*
*    BiostrucAlignAsnRead()
*
**************************************************/
NLM_EXTERN 
BiostrucAlignPtr LIBCALL
BiostrucAlignAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   BiostrucAlignPtr ptr;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* BiostrucAlign ::= (self contained) */
      atp = AsnReadId(aip, amp, BIOSTRUC_ALIGN);
   } else {
      atp = AsnLinkType(orig, BIOSTRUC_ALIGN);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = BiostrucAlignNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == BIOSTRUC_ALIGN_master) {
      ptr -> master = BiostrucAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_ALIGN_slaves) {
      ptr -> slaves = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) BiostrucAsnRead, (AsnOptFreeFunc) BiostrucFree);
      if (isError && ptr -> slaves == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_ALIGN_alignments) {
      ptr -> alignments = BiostrucAnnotSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_ALIGN_sequences) {
      ptr -> sequences = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqEntryAsnRead, (AsnOptFreeFunc) SeqEntryFree);
      if (isError && ptr -> sequences == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_ALIGN_seqalign) {
      ptr -> seqalign = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqAnnotAsnRead, (AsnOptFreeFunc) SeqAnnotFree);
      if (isError && ptr -> seqalign == NULL) {
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
   ptr = BiostrucAlignFree(ptr);
   goto ret;
}



/**************************************************
*
*    BiostrucAlignAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
BiostrucAlignAsnWrite(BiostrucAlignPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, BIOSTRUC_ALIGN);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> master != NULL) {
      if ( ! BiostrucAsnWrite(ptr -> master, aip, BIOSTRUC_ALIGN_master)) {
         goto erret;
      }
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> slaves, (AsnWriteFunc) BiostrucAsnWrite, aip, BIOSTRUC_ALIGN_slaves, BIOSTRUC_ALIGN_slaves_E);
   if (ptr -> alignments != NULL) {
      if ( ! BiostrucAnnotSetAsnWrite(ptr -> alignments, aip, BIOSTRUC_ALIGN_alignments)) {
         goto erret;
      }
   }
   AsnGenericChoiceSeqOfAsnWrite(ptr -> sequences, (AsnWriteFunc) SeqEntryAsnWrite, aip, BIOSTRUC_ALIGN_sequences, BIOSTRUC_ALIGN_sequences_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> seqalign, (AsnWriteFunc) SeqAnnotAsnWrite, aip, BIOSTRUC_ALIGN_seqalign, BIOSTRUC_ALIGN_seqalign_E);
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
*    BiostrucAlignSeqNew()
*
**************************************************/
NLM_EXTERN 
BiostrucAlignSeqPtr LIBCALL
BiostrucAlignSeqNew(void)
{
   BiostrucAlignSeqPtr ptr = MemNew((size_t) sizeof(BiostrucAlignSeq));

   return ptr;

}


/**************************************************
*
*    BiostrucAlignSeqFree()
*
**************************************************/
NLM_EXTERN 
BiostrucAlignSeqPtr LIBCALL
BiostrucAlignSeqFree(BiostrucAlignSeqPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericChoiceSeqOfFree(ptr -> sequences, (AsnOptFreeFunc) SeqEntryFree);
   AsnGenericUserSeqOfFree(ptr -> seqalign, (AsnOptFreeFunc) SeqAnnotFree);
   return MemFree(ptr);
}


/**************************************************
*
*    BiostrucAlignSeqAsnRead()
*
**************************************************/
NLM_EXTERN 
BiostrucAlignSeqPtr LIBCALL
BiostrucAlignSeqAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   BiostrucAlignSeqPtr ptr;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* BiostrucAlignSeq ::= (self contained) */
      atp = AsnReadId(aip, amp, BIOSTRUC_ALIGN_SEQ);
   } else {
      atp = AsnLinkType(orig, BIOSTRUC_ALIGN_SEQ);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = BiostrucAlignSeqNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == BIOSTRUC_ALIGN_SEQ_sequences) {
      ptr -> sequences = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqEntryAsnRead, (AsnOptFreeFunc) SeqEntryFree);
      if (isError && ptr -> sequences == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_ALIGN_SEQ_seqalign) {
      ptr -> seqalign = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqAnnotAsnRead, (AsnOptFreeFunc) SeqAnnotFree);
      if (isError && ptr -> seqalign == NULL) {
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
   ptr = BiostrucAlignSeqFree(ptr);
   goto ret;
}



/**************************************************
*
*    BiostrucAlignSeqAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
BiostrucAlignSeqAsnWrite(BiostrucAlignSeqPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, BIOSTRUC_ALIGN_SEQ);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   AsnGenericChoiceSeqOfAsnWrite(ptr -> sequences, (AsnWriteFunc) SeqEntryAsnWrite, aip, BIOSTRUC_ALIGN_SEQ_sequences, BIOSTRUC_ALIGN_SEQ_sequences_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> seqalign, (AsnWriteFunc) SeqAnnotAsnWrite, aip, BIOSTRUC_ALIGN_SEQ_seqalign, BIOSTRUC_ALIGN_SEQ_seqalign_E);
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
*    BiostrucSeqNew()
*
**************************************************/
NLM_EXTERN 
BiostrucSeqPtr LIBCALL
BiostrucSeqNew(void)
{
   BiostrucSeqPtr ptr = MemNew((size_t) sizeof(BiostrucSeq));

   return ptr;

}


/**************************************************
*
*    BiostrucSeqFree()
*
**************************************************/
NLM_EXTERN 
BiostrucSeqPtr LIBCALL
BiostrucSeqFree(BiostrucSeqPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   BiostrucFree(ptr -> structure);
   AsnGenericChoiceSeqOfFree(ptr -> sequences, (AsnOptFreeFunc) SeqEntryFree);
   return MemFree(ptr);
}


/**************************************************
*
*    BiostrucSeqAsnRead()
*
**************************************************/
NLM_EXTERN 
BiostrucSeqPtr LIBCALL
BiostrucSeqAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   BiostrucSeqPtr ptr;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* BiostrucSeq ::= (self contained) */
      atp = AsnReadId(aip, amp, BIOSTRUC_SEQ);
   } else {
      atp = AsnLinkType(orig, BIOSTRUC_SEQ);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = BiostrucSeqNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == BIOSTRUC_SEQ_structure) {
      ptr -> structure = BiostrucAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_SEQ_sequences) {
      ptr -> sequences = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqEntryAsnRead, (AsnOptFreeFunc) SeqEntryFree);
      if (isError && ptr -> sequences == NULL) {
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
   ptr = BiostrucSeqFree(ptr);
   goto ret;
}



/**************************************************
*
*    BiostrucSeqAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
BiostrucSeqAsnWrite(BiostrucSeqPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, BIOSTRUC_SEQ);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> structure != NULL) {
      if ( ! BiostrucAsnWrite(ptr -> structure, aip, BIOSTRUC_SEQ_structure)) {
         goto erret;
      }
   }
   AsnGenericChoiceSeqOfAsnWrite(ptr -> sequences, (AsnWriteFunc) SeqEntryAsnWrite, aip, BIOSTRUC_SEQ_sequences, BIOSTRUC_SEQ_sequences_E);
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
*    BiostrucSeqsNew()
*
**************************************************/
NLM_EXTERN 
BiostrucSeqsPtr LIBCALL
BiostrucSeqsNew(void)
{
   BiostrucSeqsPtr ptr = MemNew((size_t) sizeof(BiostrucSeqs));

   return ptr;

}


/**************************************************
*
*    BiostrucSeqsFree()
*
**************************************************/
NLM_EXTERN 
BiostrucSeqsPtr LIBCALL
BiostrucSeqsFree(BiostrucSeqsPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   BiostrucFree(ptr -> structure);
   AsnGenericChoiceSeqOfFree(ptr -> sequences, (AsnOptFreeFunc) SeqEntryFree);
   AsnGenericUserSeqOfFree(ptr -> seqalign, (AsnOptFreeFunc) SeqAnnotFree);
   return MemFree(ptr);
}


/**************************************************
*
*    BiostrucSeqsAsnRead()
*
**************************************************/
NLM_EXTERN 
BiostrucSeqsPtr LIBCALL
BiostrucSeqsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   BiostrucSeqsPtr ptr;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* BiostrucSeqs ::= (self contained) */
      atp = AsnReadId(aip, amp, BIOSTRUC_SEQS);
   } else {
      atp = AsnLinkType(orig, BIOSTRUC_SEQS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = BiostrucSeqsNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == BIOSTRUC_SEQS_structure) {
      ptr -> structure = BiostrucAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_SEQS_sequences) {
      ptr -> sequences = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqEntryAsnRead, (AsnOptFreeFunc) SeqEntryFree);
      if (isError && ptr -> sequences == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == BIOSTRUC_SEQS_seqalign) {
      ptr -> seqalign = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqAnnotAsnRead, (AsnOptFreeFunc) SeqAnnotFree);
      if (isError && ptr -> seqalign == NULL) {
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
   ptr = BiostrucSeqsFree(ptr);
   goto ret;
}



/**************************************************
*
*    BiostrucSeqsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
BiostrucSeqsAsnWrite(BiostrucSeqsPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objmimeAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, BIOSTRUC_SEQS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> structure != NULL) {
      if ( ! BiostrucAsnWrite(ptr -> structure, aip, BIOSTRUC_SEQS_structure)) {
         goto erret;
      }
   }
   AsnGenericChoiceSeqOfAsnWrite(ptr -> sequences, (AsnWriteFunc) SeqEntryAsnWrite, aip, BIOSTRUC_SEQS_sequences, BIOSTRUC_SEQS_sequences_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> seqalign, (AsnWriteFunc) SeqAnnotAsnWrite, aip, BIOSTRUC_SEQS_seqalign, BIOSTRUC_SEQS_seqalign_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

