/*   salutil.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  salutil.c
*
* Author:  Colombe Chappey , Hugues Sicotte
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.87 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <sequtil.h>
#include <sqnutils.h>
#include <ncbi.h>
#include <salutil.h>
#include <salstruc.h>
#include <salsap.h>
#include <edutil.h>
#include <txalign.h>
#include <bandalgn.h>
#include <blast.h>
#include <objalign.h>

#ifdef SALSA_DEBUG
#include <simutil.h>
#include <toasn3.h>
#include <utilpub.h>
#include <tfuns.h>
#endif

#define OBJ_VIRT 254

#define BAND_LIMIT 100

static SeqAlignPtr AlignmRNA2genomicToSeqAlign (SeqLocPtr slp1, SeqLocPtr slp2, SeqAlignPtr salpblast, MashPtr msp);

static Uint2 OBJ_ (Uint2 feattype)
{
  if ( feattype == FEATDEF_BAD ) return OBJ_BIOSEQ;
  return OBJ_SEQFEAT;
}

static void salutil_seqalign_write (SeqAlignPtr salp, CharPtr name)
{
  SeqAnnotPtr annot;
  AsnIoPtr aip;

  if (salp!=NULL)
  {
        annot = SeqAnnotNew();
        if (annot==NULL)
           return;
        annot->type = 2;
        annot->data = salp;

        aip = AsnIoOpen(name, "w");
        if(aip !=NULL)
        {
                        SeqAnnotAsnWrite(annot, aip, NULL);
                        AsnIoClose(aip);
        }
  }
}

/****************************************************
***   Read SeqPort.bsp from SeqPort.start to stop
***   in : spp, from + to in seq coordinates
***   out: length of buffer + buffer
***        compiled with -D SAP
*****************************************************/
extern Int4 ReadBufferFromSep (SeqPortPtr spp, CharPtr buffer, Int4 from, Int4 to, Int4 buffsegstart)
{
  Uint1    residue;
  Int4     k;
  Int4     pos;

  SeqPortSeek (spp, from, SEEK_SET);
  k = buffsegstart;
  pos = from;
  residue = SeqPortGetResidue(spp);
  while (pos < to && residue != SEQPORT_EOF)
  {
    if ( ! IS_residue(residue)) {
      /*
      switch (residue)
      {  
           case SEQPORT_VIRT:
              Message(MSG_OK,"SEQPORT_VIRT [%d=%c] at %ld\n", (int)residue, (char)residue, (long)pos);
              break;
           case SEQPORT_EOS:
              Message(MSG_OK,"[EOS]\n");
              break;
           default:
              Message(MSG_OK,"unknown char\n");
              break;
      }  
      pos++;
      */
    } else {
      buffer[k] = (Char) residue;
      k++;  
      pos++;
    }
    residue = SeqPortGetResidue(spp);
  }
  buffer[k] = '\0';
  return k;
}
 

/*****************************************************************************
***  ReadBufferFromSap 
******************************************************************************/
extern CharPtr ReadBufferFromSap (CharPtr str, CharPtr buffer, SeqAlignPtr salp, SeqIdPtr sip, Int4 from, Int4 to, Int4 *ncar)
{
  CompSegPtr  dsp;
  BoolPtr     dspstart;
  Int4Ptr     dsplens;
  Int4        dspfrom;
  Int4        sumstart;
  Int4        bufflen, buffstart;
  Int4        pre_from;
  Int4        bsplens;
  Int4        ibuff;
  Int4        istr;
  Int4        j = 0;
  Int4        k = 0;
  Int2        dim;
  Int2        numseg;
  Int2        index;
  Boolean     seen = FALSE;
  Boolean     nogap;
  Boolean     ok;

  *ncar = 0;
  if (buffer == NULL || salp == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in ReadBufferFromSap [1]\n");
         return NULL;
  }
  if ( salp->segtype != COMPSEG ) {
         ErrPostEx (SEV_ERROR, 0, 0, "fatal fail in ReadBufferFromSap"); 
         return NULL; 
  }
  if ( (dsp = (CompSegPtr) salp->segs ) == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in ReadBufferFromSap [3-2]\n");
         return NULL;
  }
  if (dsp->from == NULL || dsp->starts == NULL || dsp->lens == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in ReadBufferFromSap [3-3]\n");
         return NULL;
  }
  index = position_inIdlist (sip, dsp->ids);
  if (index == 0) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in ReadBufferFromSap [3-4]");
         return NULL;
  }
  index -= 1;
  dim = dsp->dim;
  dspfrom = *(dsp->from + index);
  dspstart = dsp->starts + index;
  dsplens = dsp->lens;
  seen = LocateInSeqAlign (from, dim, dsp->numseg, &dspstart, &dsplens, &numseg, 
                           &pre_from, &bsplens);
  if (!seen) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in ReadBufferFromSap [4]\n");
         return NULL;
  }
  ibuff = 0;
  istr = 0;
  if ( !( ok = (ibuff < to - from)) ) {
        ErrPostEx (SEV_ERROR, 0, 0, "fail in ReadBufferFromSap [6] %ld %ld %ld\n", (long) ibuff, (long) to, (long) from);
        return NULL;
  }
  nogap = (Boolean)( *dspstart );
  if ( nogap )
    buffstart= (Int4)(dspfrom +bsplens +pre_from);
  else  buffstart = 0;
  bufflen = (Int4)MIN((Int4)(*dsplens -pre_from),(Int4)(to-from +1));
  sumstart = 0;
  while ( ok ) 
  {
     if ( nogap ) 
     {
        for (j = 0; j < bufflen; j++, ibuff++, istr++) {
           buffer[ibuff] = str [istr];
           k++;
        }
        buffer[ibuff] = '\0';
        bsplens += *dsplens;
     } 
     else 
     {
        for (j = 0; j < bufflen; j++, ibuff++) {
           buffer[ibuff] = '-';
        }
        buffer[ibuff] = '\0';
     }
     numseg++;
     if (numseg > dsp->numseg) break;
     dspstart += dim; 
     dsplens++;
     nogap = (Boolean)(*dspstart);
     if ( nogap )
        buffstart = (Int4) (dspfrom +bsplens +pre_from);
     else buffstart = 0;
     bufflen = (Int4) MIN ((Int4)(*dsplens), (Int4) (to -from -ibuff +1));
  }
  *ncar = k;
  return buffer;
}


/*------------------------------------------
--  String Functions
--------------------------------------------*/
extern Boolean CCStrToInt (CharPtr str, Int2Ptr intval)
{
  Char     ch;
  Int2     i;
  Int2     len;
  Char     local [64];
  Boolean  nodigits;
  Boolean  rsult;
  int      val;
 
  rsult = FALSE;
  if (intval != NULL) {
    *intval = (Int2) 0;
  }
  len = (Int2) StringLen (str);
  if (len != 0) {
    rsult = TRUE;
    nodigits = TRUE;
    for (i = 0; i < len; i++) {
      ch = str [i];
      if (ch == ' ' || ch == '+' || ch == '-') {
      } else if (ch < '0' || ch > '9') {
        rsult = FALSE;
      } else {
        nodigits = FALSE;
      }
    }
    if (nodigits) {
      rsult = FALSE;
    }
    if (rsult && intval != NULL) {
      StringNCpy_0 (local, str, sizeof (local));
      if (sscanf (local, "%d", &val) == 1) {
        *intval = (Int2) val;
      }
    }
  }
  return rsult;
}
 
extern Boolean CCStrToLong (CharPtr str, Int4Ptr longval)
{
  Char     ch;
  Int2     i;
  Int2     len;
  Char     local [64];
  Boolean  nodigits;
  Boolean  rsult;
  long     val;
 
  rsult = FALSE;
  if (longval != NULL) {
    *longval = (Int4) 0;
  }
  len = (Int2) StringLen (str);
  if (len != 0) {
    rsult = TRUE;
    nodigits = TRUE;
    for (i = 0; i < len; i++) {
      ch = str [i];
      if (ch == ' ' || ch == '+' || ch == '-') {
      } else if (ch < '0' || ch > '9') {
        rsult = FALSE;
      } else {
        nodigits = FALSE;
      }
    }
    if (nodigits) {
      rsult = FALSE;
    }
    if (rsult && longval != NULL) {
      StringNCpy_0 (local, str, sizeof (local));
      if (sscanf (local, "%ld", &val) == 1) {
        *longval = (Int4) val;
      }
    }
  }
  return rsult;
}

extern CharPtr dashedstring (CharPtr buffer, Int4 maxbufflength)
{
  Int4 j;

  for (j = 0; j < maxbufflength; j++) 
     buffer[j] = '-';
  buffer[0] = '\0';
  return buffer;
}

extern CharPtr emptystring (CharPtr buffer, Int4 maxbufflength)
{
  Int4 j;

  if (buffer == NULL)
     return NULL;
  for (j = 0; j < maxbufflength; j++) buffer[j] = ' ';
  buffer[0] = '\0';
  return buffer;
}

extern Boolean not_empty_string (CharPtr str, Int4 lens)
{
  CharPtr  strtmp;
  Int4     j = 0;

  if (str == NULL) return FALSE;
  strtmp = str;
  for (; *strtmp != '\0' && *strtmp == ' ' && j < lens; j++, strtmp++ ) 
         continue;
  if ( j == lens || *strtmp == '\0' ) return FALSE;
  return TRUE;
}

extern Boolean stringhasnochar (CharPtr str, Int4 from, Int4 to)
{
  Char  ch;
 
  if (str != NULL) {
     if (from <= -1) 
        from = 0;
     if (to <= -1 || to > StringLen(str)) 
        to = StringLen(str);
     if (from < to && from < StringLen(str)) {
        ch = TO_UPPER(str[from]);
        while (ch != '\0' && from < to) 
        {
           if (ch >= 'A' && ch <= 'Z') {
              return FALSE;
           }
           from++;
           ch = TO_UPPER(str[from]);
        }
     }
  }
  return TRUE;
}  

extern Int1 getgapsfromstring (CharPtr str, Int4 from, Int4 to, BoolPtr *gapline)
{
  Char    ch;
  BoolPtr boolgap;
  Int4    nchar = 0; 
  Int4    len;

  if (str != NULL) {
     boolgap = *gapline;
     if (from <= -1) 
        from = 0;
     if (to <= -1 || to > StringLen(str)) 
        to = StringLen(str);
     else 
        to = MIN (to, (Int4)StringLen(str));
     len = (to-from);
     if (from < to) {
        ch = TO_UPPER(str[from]);
        while (ch != '\0' && from < to) 
        {
           if (ch >= 'A' && ch <= 'Z') {
              nchar++;
              *boolgap = FALSE;
           }
           else 
              *boolgap = TRUE;
           from++;
           ch = TO_UPPER(str[from]);
           boolgap++;
        }
        if (nchar == 0)
           return LINE_ONLYGAP;
        if (nchar < len)
           return LINE_WITHGAP;
        return LINE_NOGAP;
     }
  }
  return LINE_ONLYGAP;
}  


extern Boolean stringhasnocharplus (CharPtr str)
{
  Char  ch;
 
  if (str != NULL) {
    ch = TO_UPPER(*str);
    while (ch != '\0') {
      if ((ch >= 'A' && ch <= 'Z') || ch == '-') {
        return FALSE;
      }
      str++;
      ch = TO_UPPER(*str);
    }
  }
  return TRUE;
}  

extern CharPtr purge_string (CharPtr *st)
{
  CharPtr str;
  Int4    i, j, k, n;
  Int4    lens;

  str = *st;
  if (str==NULL)
     return NULL;
  lens = StringLen (str);
  for (i =0; i <lens; i++) {
         if (str[i] == ' ') 
         {
                j =1;
                while (str[i+j] == ' ') j++;
                lens--;
                for (k =j, n =0; k <lens; k++, n++) {
                       str[i+n] = str[k];
                }
         }
  }
  str[lens] = '\0';
  return str;
}

extern CharPtr reverse_string (CharPtr str)
{
  Char    car;
  Int4    length;
  Int4    j;

  if (str==NULL)
     return NULL;
  length = StringLen (str);
  for (j = 0; j < length / 2; j++)
  {
     car = str[j]; str[j] = str[length-j-1]; str[length-j-1] = car;
  }
  return str;
}

extern CharPtr to_lower (CharPtr str)
{
  CharPtr tmp;
 
  if (str==NULL)
     return NULL;
  tmp = str;
  while (*tmp!='\n' && *tmp!='\0') {
     *tmp = TO_LOWER(*tmp);
     tmp++;
  }
  return str;
}

/******************************************************************
****      complement_string
****                
*******************************************************************/
extern CharPtr complement_string (CharPtr str)
{
  CharPtr strp;

  for (strp = str; *strp != '\0'; strp++) {
         if (*strp == 'a') *strp = 't';
         else if (*strp == 't') *strp = 'a';
         else if (*strp == 'c') *strp = 'g';
         else if (*strp == 'g') *strp = 'c';
  }
  *strp = '\0';
  return str;
}

/*************************************************************
***
*** compare_string
***    compare 2 string by sliding in ONE direction only
***    str2 slides along str1
***    returns the position where the PERCENT of matches and
***    sequence length compared are greater
***
***************************************************************/
extern Int4 compare_string (CharPtr str1, CharPtr str2, Int4 *bestscorep)
{
  Int4  len1, len2;
  Int4  length = 0;
  Int4  score = 0;
  Int4  ratio;
  Int4  longer = 0;
  Int4  best_ratio = 0;
  Int4  best_pos = -1;
  Int4  j, k;

  if (str1 == NULL || str2 == NULL)
     return -1;
  len1 = StringLen (str1);
  len2 = StringLen (str2);
  for (j=0; j<len1; j++) {
     length = 0;
     score = 0;
     for (k=0; k<MIN(len2, len1-j); k++) {
        length ++;
        if (str1[j + k] == str2[k])
           score ++;
     }   
     ratio = (Int4) (100*score/length);
     if (ratio > best_ratio || (ratio == best_ratio && length > longer)) 
     {
        best_ratio = ratio;
        longer = length;
        best_pos = j;
     }
  }
  if (bestscorep != NULL) {
     *bestscorep = best_ratio;
  }
  return best_pos;
}  

/*************************************************************
***
*** load_seq_data
***    loads bioseq sequence into a string
***    sip: SeqId of the bioseq
***    from, to: included bondaries 
***    returns the length of the string (lenp)
***
***************************************************************/
extern CharPtr load_seq_data (SeqIdPtr sip, Int4 from, Int4 to, Boolean is_prot, Int4 *lenp)
{
  BioseqPtr        bsp;
  SeqLocPtr        slp;
  SeqPortPtr       spp;
  CharPtr          str = NULL;
  Int4             lens;

  if (from > -1 && to > -1 && from >= to)
     return NULL;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     if (from < 0 || from > bsp->length -1)
        from = 0;
     if (to < 0 || to > bsp->length -1)
        to = bsp->length -1;
     BioseqUnlock (bsp);
     if (from < to)
        slp = SeqLocIntNew (from, to, Seq_strand_plus, sip);
     else 
        slp = SeqLocIntNew (to, from, Seq_strand_minus, sip);
     if (is_prot)
        spp = SeqPortNewByLoc (slp, Seq_code_ncbieaa);
     else
        spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
     if (spp != NULL) {
        str = MemNew ((to-from+4) * sizeof(Char));
        lens = ReadBufferFromSep (spp, str, 0, to -from +1, 0);
        SeqPortFree (spp);
        if (lenp != NULL)
           *lenp = lens;
     }   
     SeqLocFree (slp);
  }
  return str;
}


extern Boolean IS_protSeqLoc (SeqLocPtr slp)
{
  CharPtr seq;
  Int4 len;

  seq = load_seq_data (SeqLocId(slp), SeqLocStart(slp), SeqLocStop(slp), TRUE, &len);
  return FALSE;
}

/*************************************************************
***
***  StringToSeqEntry :
***    in:  sequence (CharPtr) + name + length of the alignment
***    out: SeqEntryPtr
***
***************************************************************/
extern SeqEntryPtr StringToSeqEntry (CharPtr str, SeqIdPtr sip, Int4 length_align, Uint1 mol_type)
{
  SeqEntryPtr  sep;
  BioseqPtr    bsp;
  ByteStorePtr bs;
  Char         ch;

  if (str == NULL || sip == NULL) return NULL;
  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  bsp = BioseqNew ();
  if (bsp == NULL) { 
         ValNodeFree (sep); 
         return NULL; 
  }
  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;
  bsp->id = SeqIdDup (sip);
  bsp->id->next = NULL;
  SeqMgrAddToBioseqIndex (bsp);
  bsp->repr = Seq_repr_raw;
  if ( ISA_na (mol_type) ) {
         bsp->mol = Seq_mol_na;
         bsp->seq_data_type = Seq_code_iupacna;
  } 
  else {
         bsp->mol = Seq_mol_aa;
         bsp->seq_data_type = Seq_code_ncbieaa;
  }
  bsp->length = 0;
  bs = BSNew (length_align);
  bsp->seq_data = bs;
  while ((ch = *str) != '\0' && ch != ';' && bsp->length < length_align) {   
         ch = TO_UPPER (ch); 
         if ( ISA_na (mol_type) ) {
                if (ch == 'U') ch = 'T';
                if ( StringChr ("EFIJLOPQXZ-.*", ch) == NULL ) { 
                       BSPutByte ( bs, (Int2) ch );
                       bsp->length++;
                }
         } 
         else {
                if ( StringChr("JO-.", ch) == NULL ) { 
                       BSPutByte ( bs, (Int2) ch );
                       bsp->length++;
                }
         } 
         str++;
  }
  if ( bsp->length == 0 ) {
         BSFree (bs);
         BioseqFree (bsp);
         ValNodeFree (sep);
         return NULL;
  }
  return sep;
}

/*******************************************
***
***
********************************************/
extern ValNodePtr ValNodeFind (ValNodePtr head, Int2 start, Int2 index)
{
  Int2       j;

  if ( head == NULL ) return NULL;
  for (j = start; j < index && head != NULL; j++, head = head->next) 
          continue;
  return head;
}

extern ValNodePtr ValNodeFreeType (ValNodePtr *head, Int2 seqtype)
{
  ValNodePtr   vnp;

  if ( ( vnp = *head ) == NULL ) return NULL;
  for ( ; vnp != NULL; vnp = vnp->next ) {
         switch ( seqtype ) {
           case TypeEmpty: { }
               break;
           case TypeSeqInt :
               SeqIntFree ((SeqIntPtr) vnp->data.ptrvalue);
               break;
           case TypeSeqId :
               SeqIdFree ((SeqIdPtr) vnp->data.ptrvalue);
               break;
           case TypeSeqLoc :
               SeqLocFree ((SeqLocPtr) vnp->data.ptrvalue);
               break;
           case TypeSelStruct :
               SelStructDel ((SelStructPtr) vnp->data.ptrvalue);
               break;
           case TypeSelEdStruct :
               SelEdStructDel ((SelEdStructPtr) vnp->data.ptrvalue);
               break;
           default:    break;
         }
         vnp->data.ptrvalue = NULL;
  }
  ValNodeFree (vnp);
  vnp = NULL;
  return NULL;
}

extern SeqLocPtr seqloc2fuzzloc(SeqLocPtr slp, Boolean is_from, Boolean is_to)
{
   IntFuzzPtr fuzz;
   SeqIntPtr sint;

	   sint = (SeqIntPtr)slp->data.ptrvalue;
	   if(is_from){
	   	sint->if_from = IntFuzzNew();
	   	fuzz = sint->if_from;
	   	fuzz->choice = 4;
	   	fuzz->a =2;
	   }
	   if(is_to){
	     	sint->if_to = IntFuzzNew();
	     	fuzz = sint->if_to;
	     	fuzz->choice =4;
	     	fuzz->a =1;
	   }
	   return slp;
}

extern TextAlignBufPtr TextAlignBufFind (ValNodePtr anpvnp, Uint2 entityID, Uint2 itemID, Uint2 itemtype)
{
  ValNodePtr      vnp;
  TextAlignBufPtr tap;
  Uint2           tentityID;

  if ( anpvnp == NULL ) return NULL;
  for (vnp = anpvnp; vnp != NULL; vnp = vnp->next) 
  {
         tap = (TextAlignBufPtr) vnp->data.ptrvalue; 
         if ( tap != NULL)
         {
            if (OBJ_(tap->feattype) == OBJ_BIOSEQ) tentityID = tap->seqEntityID;
            else tentityID = tap->entityID;
            if (tentityID == entityID && tap->itemID == itemID && OBJ_(tap->feattype) == itemtype)
               break;
         }
  }
  if (vnp==NULL) return NULL;
  return tap;
}

extern CharPtr PNTR buf2array (ValNodePtr list, Int2 seq1, Int2 seq2)
{
  CharPtr  PNTR   tabp = NULL;
  ValNodePtr      listmp;
  TextAlignBufPtr tdp;
  Int2            nrows, j;
 
  nrows = seq2-seq1+1;
  if (nrows > 1 )  {
     tabp = (CharPtr PNTR) MemNew ((size_t)((nrows+4)*sizeof(CharPtr)));
     if (tabp != NULL)  {
        for (j=0; j<nrows+4; j++)
           tabp[j] = NULL;
        j = 0;    
        for (listmp=list; listmp!=NULL; listmp=listmp->next)
        {
           tdp = (TextAlignBufPtr) listmp->data.ptrvalue;
           if (tdp!=NULL) { 
              if (OBJ_(tdp->feattype) == OBJ_BIOSEQ) { 
                 if (tdp->buf != NULL && (j>=seq1 && j<=seq2)) {
                    tabp[j-seq1] = tdp->buf;
                    j++;
                 } 
              }  
           }  
        }
     }   
  }  
  return tabp;   
}

extern AlignNodePtr AlignNodeFind (ValNodePtr anpvnp, Uint2 entityID, Uint2 itemID, Uint2 itemtype)
{
  ValNodePtr   vnp;
  AlignNodePtr anp;

  if ( itemtype != OBJ_BIOSEQ ) return NULL;
  if ( anpvnp == NULL ) return NULL;
  for (vnp = anpvnp; vnp != NULL; vnp = vnp->next) {
         if ( (anp = (AlignNodePtr) vnp->data.ptrvalue) != NULL)
         {
              if ( anp->seq_entityID == entityID && anp->bsp_itemID == itemID )
                   break;
         }
  }
  if ( vnp == NULL ) return NULL;
  return anp;
}


extern Int2 AlignNodeIndex (ValNodePtr anpvnp, Uint2 entityID, Uint2 itemID, Uint2 itemtype)
{
  ValNodePtr   vnp;
  AlignNodePtr anp;
  Int2         index = 0;

  if ( itemtype != OBJ_BIOSEQ ) return 0;
  if ( anpvnp == NULL ) return 0;
  for (vnp = anpvnp; vnp != NULL; vnp = vnp->next, index++) {
         if ( (anp = (AlignNodePtr) vnp->data.ptrvalue) != NULL)
         {
              if ( anp->seq_entityID == entityID && anp->bsp_itemID == itemID )
                   break;
         }
  }
  if ( vnp == NULL ) return 0;
  return index;
}

/******************************************************************/
extern void OrderFeatProc (ValNodePtr vnpanp)
{
/*
  ValNodePtr   vnp;
  AlignNodePtr anp;
  AlignSegPtr  asp;

  for (vnp = vnpanp; vnp != NULL; vnp = vnp->next)
  {
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     for(asp = anp->segs; asp !=NULL; asp = asp->next)
     {
        if (asp->cnp != NULL)
           OrderCdsmRNA (&(asp->cnp));
     }
  }
*/
  return;
}

/******************************************************************/
extern ValNodePtr SeqfeatlistFree (ValNodePtr feathead)
{
  ValNodePtr     vnp, next;
  SelEdStructPtr sesp;
 
  vnp = feathead; 
  while (vnp != NULL)
  {
     next = vnp->next;
     if (vnp->data.ptrvalue != NULL) {
        sesp = (SelEdStructPtr) vnp->data.ptrvalue;
        vnp->data.ptrvalue = SelEdStructListDel (sesp);
     }
     vnp = next;
  }
  ValNodeFree (vnp);
  return NULL;
}

/******************************************************************/
extern ValNodePtr SeqfeatlistFree_fromID (ValNodePtr feathead, Uint2 entityID)
{
  ValNodePtr     vnp, pre, next;
  SelEdStructPtr sesp;
 
  vnp = feathead; 
  pre = NULL;
  while (vnp != NULL)
  {
     next = vnp->next;
     if (vnp->data.ptrvalue != NULL) {
        sesp = (SelEdStructPtr) vnp->data.ptrvalue;
        if (sesp->entityID == entityID) {
           vnp->data.ptrvalue = SelEdStructListDel (sesp);
           vnp->next = NULL;
           if (pre != NULL) pre->next = next;
           else feathead = next;
           ValNodeFree (vnp);
        }
        else pre = vnp;
     }
     vnp = next;
  }
  return feathead;
}

/******************************************************************/
extern SelEdStructPtr get_feat_fromid (ValNodePtr feat_vnp, Uint2 feattype, Uint2 ei, Uint2 ii, Int4 pos, SelEdStructPtr *prec)
{
  ValNodePtr     vnp;
  SelEdStructPtr sesp = NULL;
  SeqLocPtr      slp;
  Boolean        precedent;

  if (feat_vnp == NULL ) 
         return NULL;
  if (prec != NULL) precedent = TRUE;
  else precedent = FALSE;
  for (vnp = feat_vnp; vnp != NULL; vnp = vnp->next)
  {
         if (vnp->choice == feattype || feattype == 255) {
            sesp = (SelEdStructPtr) vnp->data.ptrvalue;
            if (sesp->entityID == ei && sesp->itemID == ii)
            {
               if (pos < 0)
                  return sesp;
               else {
                  if (precedent) *prec = NULL;
                  for (; sesp != NULL; sesp = sesp->next) {
                     if (sesp->regiontype == OM_REGION_SEQLOC && sesp->region != NULL) 
                     {
                        slp = (SeqLocPtr) sesp->region;
                        if (pos >= SeqLocStart(slp) && pos <= SeqLocStop(slp)) 
                           return sesp;
                     }
                     if (precedent) *prec = sesp;
                  }
               }
            }
         }
  }
  return NULL;
}

/******************************************************************/
extern SeqLocPtr CollectSeqLocFromAlignNode (AlignNodePtr anp)
{
  AlignSegPtr asp;
  SeqLocPtr slp = NULL;
  Int4 start, stop;
  Uint1 strand;
  Int4 current_pos;
  Int4 seglen;

  current_pos = anp->seqpos;
  if(anp->seqpos < 0)
         strand = Seq_strand_minus;
  else
         strand = Seq_strand_plus;
  for(asp = anp->segs; asp !=NULL; asp = asp->next)
  {
         seglen = 0;
         if(asp->type == REG_SEG || asp->type == DIAG_SEG)
         {
                seglen = asp->gr.right - asp->gr.left +1;
                if(strand == Seq_strand_minus)
                {
                       stop = -current_pos;
                       start = stop - (seglen-1);
                }
                else
                {
                       start = current_pos;
                       stop = start + (seglen -1);
                }
                if(slp == NULL)
                       slp = SeqLocIntNew (start, stop, strand, anp->sip);
                else
                       expand_seq_loc (start, stop, strand, slp);
         }
         current_pos += seglen;
  }
  return slp;
}

extern Int4 GetAlignLengthFromAlignNode (AlignNodePtr anp)
{
  AlignSegPtr asp;
  Int4 lens;

  lens = 0;
  for (asp = anp->segs; asp !=NULL; asp = asp->next)
  {
         if (asp->type == INS_SEG)
                lens += asp->gr.right;
         else 
                lens += asp->gr.right - asp->gr.left +1;
  }
  return lens;
}

/******************************************************************/
extern SeqIdPtr SeqIdFromAlignNode (ValNodePtr anp_lst, Uint2 entityID, Uint2 itemID, Uint2 itemtype)
{
  AlignNodePtr anp;

  anp = (AlignNodePtr) AlignNodeFind (anp_lst, entityID, itemID, itemtype);
  if (anp == NULL) 
     return NULL;
  return anp->sip;
}

extern Uint1 StrandFromAlignNode (ValNodePtr anp_lst, Uint2 entityID, Uint2 itemID, Uint2 itemtype)
{
  AlignNodePtr anp;

  anp = (AlignNodePtr) AlignNodeFind (anp_lst, entityID, itemID, itemtype);
  if (anp == NULL) 
     return Seq_strand_unknown;
  return anp->extremes.strand;
}

/*********************************************************
***
***  SeqIdPtr procedures
***    AddSeqId  : create a new seqid and add at the end
***                of the list starting with sip_head
***
***    SeqIdDupList : duplicate a list of SeqIdPtr
***
**********************************************************/
extern CharPtr matching_seqid (SeqIdPtr sip1)
{
  SeqIdPtr siptmp1, siptmp2;
  Boolean  first = TRUE;
  Char     strLog[120];
  CharPtr  str;

  for (siptmp1 = sip1; siptmp1!=NULL; siptmp1=siptmp1->next) {
    first = TRUE;
    for (siptmp2 = sip1; siptmp2!=NULL; siptmp2=siptmp2->next) {
       if (SeqIdForSameBioseq(siptmp1, siptmp2))  {
          if (first) 
             first = FALSE;
          else {
             SeqIdWrite(siptmp1, strLog, PRINTID_FASTA_LONG, 120);
             str = StringSave(strLog);
             return str;
          }
       }
    }
  }
  return NULL;
}

extern CharPtr check_seqid (Uint2 choice, CharPtr ptr)
{
  CharPtr      str;
  CharPtr      tmp;

  if (choice == SEQID_GI) {
     if (*ptr != '\0' && isdigit(*ptr)) {
        str = StringSave (ptr);
        tmp = ptr;
        tmp = StringMove (tmp, "gi|"); 
        tmp = StringMove (tmp, str); 
        *tmp = '\0';
        MemFree (str);
     }
  }
  return ptr;
}

/*********************************************************
***
***    AddSeqId  : create a new seqid and add at the end
***                of the list starting with sip_head
***
**********************************************************/
extern SeqIdPtr AddSeqId (SeqIdPtr *sip_head, SeqIdPtr sip)
{
  SeqIdPtr sip_tmp, 
           sip_copy;

  sip_copy = SeqIdDup (sip);
  sip_tmp = sip_copy->next;
  sip_copy->next = NULL;
  if (sip_tmp!=NULL)
     SeqIdFree (sip_tmp);
  if ( (sip_tmp = *sip_head) != NULL ) {
     while (sip_tmp->next != NULL) 
        sip_tmp = sip_tmp->next; 
     sip_tmp->next = sip_copy;
  }  
  else {
     *sip_head = sip_copy;
  }
  return (*sip_head);

}

/*******************************************************
***
***   SeqIdDupList : duplicate a list of SeqIdPtr
***
*******************************************************/
extern SeqIdPtr SeqIdDupList (SeqIdPtr id_list)
{
  SeqIdPtr     sip=NULL;
  SeqIdPtr     sid;

  for (sid = id_list; sid != NULL; sid = sid->next) {
         sip = AddSeqId (&sip, sid);  
  }
  return sip;
}

extern SeqIdPtr SeqIdDupBestList (SeqIdPtr id_list)
{
  SeqIdPtr     sip=NULL;
  SeqIdPtr     sid, sid2;
  BioseqPtr    bsp;

  for (sid = id_list; sid != NULL; sid = sid->next) {
     sid2 = NULL;
     bsp = BioseqLockById (sid);
     if (bsp!=NULL) {
        sid2 = SeqIdFindBest(bsp->id, 0);
        BioseqUnlock (bsp);
     }
     if (sid2!=NULL)
        sip = AddSeqId (&sip, sid2);
     else 
        sip = AddSeqId (&sip, sid);  
  }
  return sip;
}

extern SeqIdPtr SeqIdListfromSeqLoc (ValNodePtr vnpslp)
{
  SeqIdPtr     sip=NULL, siptmp;
  ValNodePtr   vnp=NULL;
  Int2         j = 0, k;
  for (vnp = vnpslp; vnp != NULL; vnp = vnp->next)
  {
         sip = AddSeqId (&sip, SeqLocId ((SeqLocPtr) vnp->data.ptrvalue));  
         j++;
  }
  if (sip!=NULL) {
     for (siptmp=sip, k=0; k<j-1; siptmp=siptmp->next, k++) continue;
     siptmp->next = NULL;
  }
  return sip;
}

extern SeqIdPtr ValNodeSeqIdListDup (ValNodePtr id_list)
{
  ValNodePtr   vnp=NULL;
  SeqIdPtr     sip = NULL;

  for (vnp = id_list; vnp != NULL; vnp = vnp->next) 
  {
         sip = AddSeqId (&sip, (SeqIdPtr) vnp->data.ptrvalue);  
  }
  return sip;
}
/*******************************************************
***
***   SeqIdListToCharArray
***
*******************************************************/
extern CharPtr PNTR SeqIdListToCharArray (SeqIdPtr id_list, Int2 n)
{
  CharPtr PNTR  idarray; 
  SeqIdPtr      idtmp; 
  Int2          k;

  idarray = MemNew ((size_t) ((n+1) * sizeof(CharPtr)));
  idtmp = id_list;
  for (k = 0; k < n && idtmp!=NULL; k++, idtmp =idtmp->next) {
         idarray [k] = (CharPtr) MemNew ((size_t) ((52) * sizeof(Char)));
         SeqIdWrite (idtmp, idarray [k], PRINTID_FASTA_SHORT, 50);
  }
  return idarray;
}

/******************************************************************/
extern Int2 position_inIdlist (SeqIdPtr a, SeqIdPtr b)
{
  SeqIdPtr siptmp;
  Int2     index = 1;
  Boolean  retval;
 
  if (a!=NULL && b!=NULL) 
     for (siptmp = b; siptmp != NULL; siptmp = siptmp->next)
     {
        retval = SeqIdForSameBioseq (a, siptmp);
        if (retval == TRUE)
           return index;
        index++;
     }
  return 0;
}

/******************************************************************/
extern SeqIdPtr SeqIdReplaceID (SeqIdPtr head, SeqIdPtr pre, SeqIdPtr sip, SeqIdPtr next)
{
  SeqIdPtr tmp;

  if (pre == NULL)
  {
     head = SeqIdDup(sip);
     head->next = next;
     return head;
  }
  tmp = pre->next;
  pre->next = NULL;
  tmp->next = NULL;
  SeqIdFree (tmp);
  pre->next = SeqIdDup(sip);
  pre->next->next = next;
  return head;
}

typedef struct ccid {
  SeqIdPtr    sip;
  SeqEntryPtr sep;
  BioseqPtr   bsp;
} CcId, PNTR CcIdPtr;

static void SeqAnnotReplaceID (SeqAnnotPtr sap, SeqIdPtr newsip)
{
  SeqFeatPtr  sfp;

  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->location != NULL)
           sfp->location = SeqLocReplaceID (sfp->location, newsip);
        if (sfp->product != NULL)
           sfp->product = SeqLocReplaceID (sfp->product, newsip);
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

extern BioseqPtr BioseqReplaceID (BioseqPtr bsp, SeqIdPtr newsip)
{
  if (bsp!=NULL)
  {
     SeqIdFree (bsp->id);
     bsp->id = newsip;
     SeqMgrReplaceInBioseqIndex (bsp);
     SeqAnnotReplaceID (bsp->annot, newsip);
  }
  return bsp;
}

static void SeqEntryReplaceSeqIDCallBack (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  SeqIdPtr PNTR      sipp;
  SeqIdPtr           newsip;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;

  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     sipp = mydata;
     newsip = *sipp;
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (ISA_na (bsp->mol)) {
           BioseqReplaceID (bsp, newsip);
        }
     } 
     else if (IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp->annot) {
           SeqAnnotReplaceID (bssp->annot, newsip);
        }
     }
  }
}
  

static void FindSeqEntryForSeqIdCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  CcIdPtr            cip;

  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     cip = (CcIdPtr)mydata;
     if (cip->sep==NULL && IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL && ISA_na (bsp->mol)) 
        {
           if (SeqIdForSameBioseq(cip->sip, bsp->id)) {
              cip->sep = sep;
              cip->bsp = bsp;
           }
        }
     }
  }
}

extern SeqEntryPtr SeqEntryReplaceSeqID (SeqEntryPtr source_sep, SeqIdPtr sip)
{
  SeqEntryPtr sep;
  SeqLocPtr   slp;
  SeqIdPtr    newsip;
  CcId        ci;
  Uint2       entityID;

  ci.sip = SeqIdDup (sip);
  ci.sep = NULL;
  ci.bsp = NULL;
  SeqEntryExplore(source_sep,(Pointer)&ci, FindSeqEntryForSeqIdCallback);
  if (ci.sep && ci.bsp)
  {
     entityID = ObjMgrGetEntityIDForPointer (ci.bsp);
     sep = GetBestTopParentForData (entityID, ci.bsp);
     slp = SeqLocIntNew (0, ci.bsp->length-1, Seq_strand_plus, ci.bsp->id);
     newsip = MakeNewProteinSeqId(slp, NULL);
     ValNodeFree (slp);
     SeqEntryExplore (sep, &newsip, SeqEntryReplaceSeqIDCallBack);
  }
  SeqIdFree (ci.sip);
  return source_sep;
}

/*********************************************************
***
***  ScorePtr procedures
***
**********************************************************/
extern ScorePtr ScoreDup (ScorePtr sp)

{
  ScorePtr sp_copy;
  sp_copy = ScoreNew();
  sp_copy->id = (ObjectIdPtr) ObjectIdDup (sp->id);
  sp_copy->choice = sp->choice;
  sp_copy->value = sp->value;
  sp_copy->next = NULL;
  return sp_copy;
}

extern ScorePtr ScoreDupAdd (ScorePtr *sp_head, ScorePtr sp)
{
  ScorePtr sp_tmp, sp_copy;

  sp_copy = ScoreDup (sp);
  if ( (sp_tmp = *sp_head) != NULL ) {
         while (sp_tmp->next != NULL) sp_tmp = sp_tmp->next; 
         sp_tmp->next = sp_copy;
  }  
  else *sp_head = sp_copy;
  return *sp_head;
}

extern ScorePtr ScoreAdd (ScorePtr *sp_head, ScorePtr sp)
{
  ScorePtr sp_tmp;

  if ( (sp_tmp = *sp_head) != NULL ) {
         while (sp_tmp->next != NULL) sp_tmp = sp_tmp->next; 
         sp_tmp->next = sp;
  }  
  else *sp_head = sp;
  return *sp_head;
}


/*********************************************************
***
***  SeqLocPtr procedures
***
**********************************************************/
extern Int2 chkloc (SeqIdPtr sip, Int4 position, ValNodePtr vnp, Int4 *newpos)
{
  ValNodePtr vnptmp;
  SeqIdPtr   siptmp;
  SeqLocPtr  slp;

  *newpos = 0;
  for (vnptmp=vnp; vnptmp!=NULL; vnptmp=vnptmp->next)
  {
     slp = (SeqLocPtr)vnptmp->data.ptrvalue;
     siptmp = SeqLocId (slp);
     if (siptmp!=NULL)
        if (SeqIdForSameBioseq (sip, siptmp)) {
           if (position >= SeqLocStart(slp) && position <= SeqLocStop(slp)) {
              *newpos = position;
              return 0;
           }
           if (position==APPEND_RESIDUE || position>=SeqLocStart(slp)+SeqLocLen(slp)) {
              *newpos = SeqLocStart(slp) + SeqLocLen(slp); 
              return APPEND_RESIDUE;
           }
           if (position < SeqLocStart(slp)) {
              *newpos = SeqLocStart(slp);
              return GAP_RESIDUE; 
           }
        }
  }
  return NO_RESIDUE;
}

extern SeqLocPtr expand_seq_loc(Int4 start, Int4 stop, Uint1 strand, SeqLocPtr loc)
{
   SeqIntPtr sint;
   SeqPntPtr spp;
 
        if(loc->choice == SEQLOC_INT)
        {
                sint = loc->data.ptrvalue;
                if(start != -1 && start < sint->from)
                        sint->from = start;
                if(stop != -1 && stop > sint->to)
                        sint->to = stop;
                if(strand != 0 && sint->strand != strand)
                        sint->strand = strand;
                loc->data.ptrvalue = sint;
        }
        else if(loc->choice == SEQLOC_PNT)
        {
                spp = (SeqPntPtr)(loc->data.ptrvalue);
                spp->point = start;
                spp->strand = strand;
                loc->data.ptrvalue = spp;
        }
 
        return loc;
}

extern Int4 MaxLengthSeqLoc (ValNodePtr sqloc_list)
{
  ValNodePtr       vnp;
  SeqLocPtr        slp;
  Int4             len, maxlen = 0;

  for (vnp = sqloc_list; vnp != NULL; vnp = vnp->next)
  { 
         slp = (SeqLocPtr) vnp->data.ptrvalue;
         if ( ( len = SeqLocLen (slp)) > maxlen ) 
                maxlen = len;
  }
  return maxlen;
}

extern Boolean SeqLocListMatch (ValNodePtr vnp1, ValNodePtr vnp2, Boolean *Fp, Boolean *Tp)
{
  ValNodePtr tmp1, tmp2;
  SeqLocPtr  slp1, slp2;
  Boolean    p5short = FALSE, p3short=FALSE;

  for (tmp1=vnp1, tmp2=vnp2; tmp1!=NULL && tmp2!=NULL; tmp1=tmp1->next, tmp2=tmp2->next)
  {
     slp1 = (SeqLocPtr) tmp1->data.ptrvalue;
     slp2 = (SeqLocPtr) tmp2->data.ptrvalue;
     if (!p5short) 
            p5short = ( SeqLocStart(slp1) < SeqLocStart(slp2) );
     if (!p3short) 
            p3short  = ( SeqLocStop (slp1) > SeqLocStop (slp2) );
     if (p5short && p3short) break;
  }
  *Fp = p5short;
  *Tp = p3short;
  if (p5short || p3short) return FALSE;
  return TRUE;
} 


/***********************************************************************
***    
***    SeqEntryToSeqLoc
***      read SeqEntry (1->Bioseq or 2->BioseqSetPtr)
***      return list of ValNodePtr->SeqLocPtr
***      The strand of the seqloc is 0, whatever the strand of the bsp 
************************************************************************/
typedef struct ccid3
{
  ValNodePtr vnp;
  Uint1      bsp_mol;
} CcId3, PNTR CcId3Ptr;

static void ListSeqEntryCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  CcId3Ptr           ccp;

  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     ccp = (CcId3Ptr)mydata;
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
         if (ccp->bsp_mol==0 || ISA_aa(ccp->bsp_mol)==ISA_aa(bsp->mol))
         {
           sip = SeqIdFindBest(bsp->id, 0);
           if (sip!=NULL) {
              slp = SeqLocIntNew (0, bsp->length - 1, Seq_strand_plus, sip);
              ValNodeAddPointer (&(ccp->vnp), 0, slp);
           }
         }
        }
     }
  }
}

extern ValNodePtr SeqEntryToSeqLoc (SeqEntryPtr sep, Int2 *n, Uint1 bsp_mol)
{
  CcId3              cc;
  ValNodePtr         vnp = NULL;
  SeqAlignPtr        salp = NULL;
  Int2               j=0;

  cc.vnp=NULL;
  cc.bsp_mol=bsp_mol;
  SeqEntryExplore(sep,(Pointer)&cc, ListSeqEntryCallback);
  for (vnp = cc.vnp;vnp!=NULL; vnp=vnp->next)
  {
     j++;
  } 
  *n = j;
  return cc.vnp;
}


/************************************************************/
extern SelStructPtr SelStructNew (Uint2 ssp_ed, Uint2 ssp_id, Uint2 ssp_it, Int4 from, Int4 to, SeqIdPtr sip, Uint1 strand, Boolean is_fuzz)
{
  SelStructPtr   ssp;
  SeqLocPtr      slp;

  ssp = (SelStructPtr) MemNew (sizeof (SelStruct));
  if (ssp == NULL) return NULL;
  ssp->entityID = ssp_ed;
  ssp->itemID = ssp_id;
  ssp->itemtype = ssp_it;
  ssp->regiontype = 0;
  ssp->region = NULL;
  if (from >= 0 && sip != NULL) {
     if (is_fuzz)
         slp = fuzz_loc (from, to, strand, sip, TRUE, TRUE);
     else 
         slp = SeqLocIntNew (from, to, strand, sip);
     if ( slp != NULL ) {
         ssp->regiontype = OM_REGION_SEQLOC;
         ssp->region = (Pointer) slp;
     }
  }
  ssp->prev = NULL;
  ssp->next = NULL;
  return ssp;
}
/******************************************************************/
extern SelStructPtr SelStructCpy (SelStructPtr ssp, SelStructPtr ssp2)
{
  SeqLocPtr  slp, slp2;
  if (ssp == NULL) return NULL;
  ssp2->entityID = ssp->entityID;
  ssp2->itemID = ssp->itemID;
  ssp2->itemtype = ssp->itemtype;
  if ( ssp->regiontype == OM_REGION_SEQLOC) {
         ssp2->regiontype = OM_REGION_SEQLOC;
         if ( ssp->region == NULL ) {
               ErrPostEx (SEV_WARNING, 0, 0, "Fail in SelStructCpy [1]");
               return NULL;
         }
         if ( ssp2->region != NULL ) 
               SeqLocFree ((SeqLocPtr)ssp2->region);
         slp = (SeqLocPtr) ssp->region;
         slp2 = SeqLocIntNew (SeqLocStart(slp), SeqLocStop(slp), 
                              SeqLocStrand(slp), SeqLocId(slp));
         if ( slp2 == NULL ) {
               ErrPostEx(SEV_WARNING, 0, 0, "Fail in SelStructCpy [2]");
               return NULL;
         }
         ssp2->region = (Pointer) slp2;
  } 
  else {
         ssp2->regiontype = 0;
         ssp2->region = NULL;
  }
  ssp2->next = NULL;
  ssp2->prev = NULL;
  return ssp2;
}

/******************************************************************/
extern SelStructPtr SelStructDup (SelStructPtr ssp)
{
  SelStructPtr ssp2;

  if (ssp == NULL) return NULL;
  ssp2 = (SelStructPtr) MemNew (sizeof (SelStruct));
  ssp2 = SelStructCpy (ssp, ssp2);
  return ssp2;
}

/******************************************************************/
extern SelStructPtr SelStructAdd (SelStructPtr head, SelStructPtr ssp)
{
  SelStructPtr ssptmp;
  if (head == NULL) return ssp;
  for (ssptmp = head; ssptmp->next != NULL; ssptmp=ssptmp->next) continue;
  ssptmp->next = ssp;
  ssp->prev = ssptmp;
  return head;
}

/******************************************************************/
extern SelStructPtr SelStructDel (SelStructPtr ssp)
{

  if (ssp == NULL) return NULL;
  ssp->next = NULL;
  if ( ssp->region != NULL ) 
         SeqLocFree ((SeqLocPtr) ssp->region);
  MemFree (ssp);
  return NULL;
}

/******************************************************************/
extern SelStructPtr SelStructDelList (SelStructPtr ssp)
{
  SelStructPtr tmp, next;
 
  if (ssp!=NULL) {
     tmp=ssp;
     while (tmp!=NULL) {
        next = tmp->next;
        tmp->next = NULL;
        SelStructDel (tmp);
        tmp=next;
     }
  }   
  return NULL;
}

/*****************************************************************/
extern void setposition_tossp (SelStructPtr ssp, Int4 from, Int4 to)
{
  SeqLocPtr        slp;
  SeqIntPtr        sitcaret;
  if (ssp == NULL) return;
  if (ssp->region == NULL) return;
  slp = (SeqLocPtr) ssp->region;
  sitcaret = (SeqIntPtr) slp->data.ptrvalue;
  sitcaret->from = from;
  sitcaret->to = to;
}

/******************************************************************/
extern Boolean is_samessp (SelStructPtr ssp1, SelStructPtr ssp2)
{
  if (ssp1 == NULL || ssp2 == NULL) return FALSE;
  if (ssp1->entityID != ssp2->entityID) return FALSE;
  if (ssp1->itemID != ssp2->itemID) return FALSE;
  if (ssp1->itemtype != ssp2->itemtype) return FALSE;
  return TRUE;
}

/******************************************************************/
extern Boolean is_sameId (Uint2 sei, Uint2 sii, Uint2 sit, Uint2 sist, Uint2 ei, Uint2 ii, Uint2 it, Uint2 ist)
{
  if ( ei < 255 && sei != ei) return FALSE;
  if ( ii < 255 && sii != ii) return FALSE;
  if ( it < 255 && sit != it) return FALSE;
  if ( ist< 255 && sist!= ist)return FALSE;
  return TRUE;
}

/******************************************************************/
extern Boolean is_samepos (SelStructPtr ssp1, SelStructPtr ssp2)
{
  SeqLocPtr       slp1, slp2;

  if (ssp1->regiontype == 0 || ssp2->regiontype == 0) return FALSE;
  if (ssp1->region == NULL || ssp2->region == NULL) return FALSE;
  if ( !is_samessp (ssp1, ssp2) ) return FALSE;
  slp1 = (SeqLocPtr) ssp1->region;
  slp2 = (SeqLocPtr) ssp2->region;
  if ( slp1 == NULL || slp2 == NULL) return FALSE;
  if ( SeqLocStart(slp1) != SeqLocStart(slp2) ) return FALSE;
  if ( SeqLocStop(slp1) != SeqLocStop(slp2) ) return FALSE;
  return TRUE;
}

/******************************************************************/
extern ValNodePtr del_ssp_fromid (ValNodePtr headp, Uint2 itemsubtype, SelEdStructPtr target)
{
  ValNodePtr     vnphead, pre_vnp = NULL;
  SelEdStructPtr ssp,  pressp = NULL;
  Boolean        found_ssp = FALSE;

  if (headp == NULL || target == NULL) {
         ErrPostEx(SEV_WARNING, 0, 0, "fail in del_ssp_fromid [1]");
         return NULL;
  }
  pre_vnp = NULL;
  for  (vnphead = headp; vnphead != NULL; vnphead = vnphead->next) 
  {
     if ( itemsubtype == 255 || vnphead->choice == itemsubtype ) {
        ssp = (SelEdStructPtr) vnphead->data.ptrvalue; 
        if (is_sameses (ssp, target)) {
           pressp = NULL;
           for (; ssp != NULL; ssp = ssp->next) 
           {
              if (include_ssp ((SeqLocPtr) ssp->region, (SeqLocPtr) target->region)) 
              {
                 found_ssp = TRUE; 
                 break;
              }
              pressp = ssp;
           }
           if (found_ssp) break;
        }
     }
     pre_vnp = vnphead;
  } 
  if (vnphead == NULL || ssp == NULL) return headp;
  SelEdStructListDel ((SelEdStructPtr) vnphead->data.ptrvalue);
  vnphead->data.ptrvalue = NULL;
  if (pre_vnp == NULL) {
     if (vnphead->next == NULL) {
        headp = NULL;
     }
     else {
        headp = vnphead->next;
        vnphead->next = NULL;
     }
  }
  else if (vnphead->next == NULL) {
     pre_vnp->next = NULL;
  }
  else {
     pre_vnp->next = vnphead->next;
     vnphead->next = NULL;
  }
  ValNodeFree (vnphead);
  return headp;
}

/******************************************************************/
extern Boolean include_ssp (SeqLocPtr slp1, SeqLocPtr slp2)
{
  if ( SeqLocStart(slp1) <= SeqLocStart(slp2) 
  && SeqLocStop(slp1) >= SeqLocStop(slp2) ) return TRUE;
  return FALSE;
}

/******************************************************************/
extern Int4 overlapp_startssp (SeqLocPtr slp1, SeqLocPtr slp2)
{
  if ( SeqLocStart(slp1) < SeqLocStart(slp2) 
  &&  SeqLocStop(slp1) < SeqLocStop(slp2) ) 
     return (SeqLocLen(slp1) - (SeqLocStart(slp2) - SeqLocStart(slp1)));
  return 0;
}

/******************************************************************/
extern Boolean overlapp_ssp (SeqLocPtr slp1, SeqLocPtr slp2)
{
  if ( SeqLocStop(slp1) < SeqLocStart(slp2) ) return FALSE;
  if ( SeqLocStart(slp1) > SeqLocStop(slp2) ) return FALSE;
  return TRUE;
}

/******************************************************************/
extern Boolean precede_ssp (SeqLocPtr slp1, SeqLocPtr slp2)
{
  if ( SeqLocStart(slp2) >= SeqLocStop(slp1) ) return TRUE;
  return FALSE;
}

/******************************************************************/
extern Boolean succeed_ssp (SeqLocPtr slp1, SeqLocPtr slp2)
{
  if ( SeqLocStop(slp2) == SeqLocStart(slp1) ) return TRUE;
  return FALSE;
}

extern SelStructPtr addssp (SelStructPtr *ssphead, Uint2 choice, Pointer pt, Uint2 iID)
{
  SelStructPtr ssp, hssp;

  ssp = (SelStructPtr) MemNew (sizeof (SelStruct));
  if (ssp == NULL) 
     return *ssphead;
  ssp->entityID = 0;
  ssp->itemID = iID;
  ssp->itemtype = choice;
  ssp->regiontype = OM_REGION_SEQLOC;  
  ssp->region = (Pointer) pt;
  ssp->next = NULL;
  ssp->prev = NULL;
  hssp = *ssphead; 
  if (hssp == NULL) {
     *ssphead = ssp;
  }
  else {
     for (; hssp->next != NULL; hssp = hssp->next) continue;
     hssp->next = ssp;
     ssp->prev = hssp;
  }
  return ssp;
}


/******************************************************************/
extern SelEdStructPtr new_seledstruct (Uint2 ssp_ed, Uint2 ssp_id, Uint2 ssp_it, Uint2 bspiID, Int4 from, Int4 to, SeqIdPtr sip, Uint1 strand, Boolean is_fuzz, CharPtr label, Pointer data, Int4 offset, Uint1 codonstart)
{
  SelEdStructPtr sesp;
  SeqLocPtr      slp;
  ValNodePtr     datavnp;

  sesp = (SelEdStructPtr) MemNew (sizeof (SelEdStruct));
  if (sesp == NULL) return NULL;
  sesp->entityID = ssp_ed;
  sesp->itemID = ssp_id;
  sesp->itemtype = ssp_it;
  sesp->bsp_itemID = bspiID;
  if (sesp->regiontype == OM_REGION_SEQLOC) 
         sesp->regiontype = 0;
  if (from < 0 || sip == NULL) {
         sesp->regiontype = 0;
         sesp->region = NULL;
  }
  else {
     sesp->region = NULL; 
     if (is_fuzz)
         slp = fuzz_loc (from, to, strand, sip, TRUE, TRUE);
     else 
         slp = SeqLocIntNew (from, to, strand, sip);
     if ( slp != NULL ) {
         sesp->regiontype = OM_REGION_SEQLOC;
         sesp->region = (Pointer) slp;
     }
     else   sesp->regiontype = 0; 
  }
  if (label != NULL && StringLen (label) > 0)
       StringCpy(sesp->label, label);
  else sesp->label[0] = '\0';
  if (data != NULL)
  {
     datavnp = ValNodeNew (NULL);
     datavnp->choice = 0;
     datavnp->data.ptrvalue = data;
     sesp->data = datavnp;
  }
  else sesp->data = NULL;
  sesp->codonstart = codonstart;
  sesp->offset = offset;
  sesp->dirty = TRUE;
  sesp->visible = TRUE;
  sesp->prev = NULL;
  sesp->next = NULL;
  return sesp;
}

extern SelEdStructPtr new_seledstruct_fromseqloc (Uint2 ssp_ed, Uint2 ssp_id, Uint2 ssp_it, Uint2 bspiID, SeqLocPtr slp, CharPtr label, Pointer data, Int4 offset, Uint1 codonstart)
{
  SelEdStructPtr sesp;
  ValNodePtr     datavnp;

  sesp = (SelEdStructPtr) MemNew (sizeof (SelEdStruct));
  if (sesp == NULL) return NULL;
  sesp->entityID = ssp_ed;
  sesp->itemID = ssp_id;
  sesp->itemtype = ssp_it;
  sesp->bsp_itemID = bspiID;
  if (sesp->regiontype == OM_REGION_SEQLOC) 
     sesp->regiontype = 0;
  if (slp == NULL) {
     sesp->regiontype = 0;
     sesp->region = NULL;
  }
  else {
     sesp->regiontype = OM_REGION_SEQLOC;
     sesp->region = (Pointer) slp;
  }
  if (label != NULL && StringLen (label) > 0)
     StringCpy(sesp->label, label);
  else sesp->label[0] = '\0';
  if (data != NULL)
  {
     datavnp = ValNodeNew (NULL);
     datavnp->choice = 0;
     datavnp->data.ptrvalue = data;
     sesp->data = datavnp;
  }
  else sesp->data = NULL;
  sesp->codonstart = codonstart;
  sesp->offset = offset;
  sesp->dirty = TRUE;
  sesp->visible = TRUE;
  sesp->prev = NULL;
  sesp->next = NULL;
  return sesp;
}

/***************************************************************
*** sesp_to_slp
***    make seqloc from  ses->region (SelEdStructPtr) 
***    after changing the alignment coordinates into 
***    sequence coordinates
***    
***    the seqloc is NOT PARTIAL
***        fuzz_loc (start, stop, strand, sip, TRUE, TRUE);
***    but COMPLETE: 
***        SeqLocIntNew (start, stop, strand, sip);)
***************************************************************/
extern SeqLocPtr sesp_to_slp (SelEdStructPtr ses, SeqAlignPtr salp, ValNodePtr sqlocs, Boolean partial)
{
  SelEdStructPtr sesp;
  SeqLocPtr      tmp, new, 
                 slp, slp1;
  SeqIdPtr       sip;
  Int4           start, stop;
  Uint1          strand;

  Int2           j,k;

  if (ses->next == NULL)
  {
     tmp = (SeqLocPtr) ses->region;
     sip = SeqLocId (tmp);
     start= (Int4)AlignCoordToSeqCoord (SeqLocStart(tmp), sip, salp, sqlocs,0);
     stop = (Int4)AlignCoordToSeqCoord (SeqLocStop(tmp), sip, salp, sqlocs, 0);
     slp = SeqLocIntNew (start, stop, SeqLocStrand (tmp), sip);
     return slp;
  }
  slp1 = (SeqLocPtr) ValNodeNew (NULL);
  slp1->choice = SEQLOC_PACKED_INT;
  sesp = ses;
  tmp = (SeqLocPtr) sesp->region;
  if (SeqLocStrand (tmp) == Seq_strand_minus)
  {
     j=0;
     while (sesp->next!=NULL) {
        j++; sesp = sesp->next;
     }
     tmp = (SeqLocPtr) sesp->region;
  }
  strand = SeqLocStrand (tmp);
  sip = SeqLocId (tmp);
  start= (Int4) AlignCoordToSeqCoord (SeqLocStart (tmp), sip, salp, sqlocs, 0);
  stop = (Int4) AlignCoordToSeqCoord (SeqLocStop (tmp), sip, salp, sqlocs, 0);
  if (partial)
     new = fuzz_loc (start, stop, strand, sip, TRUE, TRUE);
  else
     new = SeqLocIntNew (start, stop, strand, sip);
  slp1->data.ptrvalue = new;
  slp = new;
  if (strand != Seq_strand_minus)
  {
     sesp = sesp->next;
     for (; sesp != NULL; sesp = sesp->next)
     {
        tmp = (SeqLocPtr) sesp->region;
        start=(Int4)AlignCoordToSeqCoord(SeqLocStart(tmp),sip, salp, sqlocs, 0);
        stop =(Int4)AlignCoordToSeqCoord(SeqLocStop(tmp), sip, salp, sqlocs,0);
        new = SeqLocIntNew (start, stop, strand, sip);
        slp->next = new;
        slp = slp->next;
     }
  }  
  else {
     while (j>0) {
        sesp=ses;
        for (k=1; k<j; k++) sesp=sesp->next;
        tmp = (SeqLocPtr) sesp->region;
        start=(Int4)AlignCoordToSeqCoord(SeqLocStart(tmp), sip, salp, sqlocs,0);
        stop =(Int4)AlignCoordToSeqCoord(SeqLocStop(tmp), sip, salp, sqlocs, 0);
        new = SeqLocIntNew (start, stop, strand, sip);
        slp->next = new;
        slp = slp->next;
        j--;
     }
  }
  return slp1;
}


static SelEdStructPtr SelEdStructCpy (SelEdStructPtr ssp, SelEdStructPtr ssp2)
{
  SeqLocPtr  slp, slp2;
  if (ssp == NULL) return NULL;
  ssp2->entityID = ssp->entityID;
  ssp2->itemID = ssp->itemID;
  ssp2->itemtype = ssp->itemtype;
  if ( ssp->regiontype == OM_REGION_SEQLOC && ssp->region != NULL) {
         ssp2->regiontype = OM_REGION_SEQLOC;
         if ( ssp2->region != NULL ) 
               SeqLocFree ((SeqLocPtr)ssp2->region);
         slp = (SeqLocPtr) ssp->region;
         slp2 = SeqLocIntNew (SeqLocStart(slp), SeqLocStop(slp), SeqLocStrand(slp), SeqLocId(slp));
         if ( slp2 == NULL ) {
               ErrPostEx(SEV_WARNING, 0, 0, "Fail in SelStructCpy [2]");
               return NULL;
         }
         ssp2->region = (Pointer) slp2;
  } 
  else {
         ssp2->regiontype = 0;
         ssp2->region = NULL;
  }
  ssp2->data = NULL;
  ssp2->next = NULL;
  ssp2->prev = NULL;
  return ssp2;
}

/******************************************************************/
extern SelEdStructPtr SelEdStructDup (SelEdStructPtr ssp)
{
  SelEdStructPtr ssp2;

  if (ssp == NULL) return NULL;
  ssp2 = (SelEdStructPtr) MemNew (sizeof (SelEdStruct));
  ssp2 = SelEdStructCpy (ssp, ssp2);
  return ssp2;
}

/******************************************************************/
extern SelEdStructPtr SelEdStructAdd (SelEdStructPtr head, SelEdStructPtr ssp)
{
  SelEdStructPtr ssptmp;
  if (head == NULL) return ssp;
  for (ssptmp = head; ssptmp->next != NULL; ssptmp=ssptmp->next) continue;
  ssptmp->next = ssp;
  ssp->prev = ssptmp;
  return head;
}

/******************************************************************/
extern SelEdStructPtr SelEdStructDel (SelEdStructPtr ssp)
{
  ValNodePtr  vnpdata;

  if (ssp == NULL) return NULL;
  if ( ssp->data != NULL ) 
  {
           vnpdata = ssp->data;
           if ( vnpdata->data.ptrvalue != NULL ) {
              if (ssp->prev == NULL)
                 vnpdata->data.ptrvalue = MemFree(vnpdata->data.ptrvalue);
              else 
                 vnpdata->data.ptrvalue = NULL;
           }
           vnpdata = ValNodeFree (vnpdata);
  }
  if ( ssp->region != NULL ) 
         ssp->region = SeqLocFree ((SeqLocPtr) ssp->region);
  MemFree (ssp);
  return NULL;
}

/******************************************************************/
extern SelEdStructPtr SelEdStructListDel (SelEdStructPtr ssp)
{
  SelEdStructPtr next;
  ValNodePtr     vnpdata;
  Boolean        first = TRUE;

  if (ssp == NULL) return NULL;
  while (ssp != NULL) 
  {
     next = ssp->next;
     if ( ssp->region != NULL ) 
           ssp->region = SeqLocFree ((SeqLocPtr) ssp->region);
     if ( ssp->data != NULL ) 
     {
           vnpdata = ssp->data;
           if ( vnpdata->data.ptrvalue != NULL ) {
              if (first) 
		  vnpdata->data.ptrvalue = MemFree (vnpdata->data.ptrvalue);
	      else 
		  vnpdata->data.ptrvalue = NULL;
           }
           vnpdata = ValNodeFree (vnpdata);
     }
     MemFree (ssp);
     ssp = next;
	 if (first) first = FALSE;
  }
  return NULL;
}

extern void set_seqnot_visible (Uint2 eID, Uint2 iID, SelEdStructPtr sesp)
{
  SelEdStructPtr tmp;
  for (tmp=sesp;tmp!=NULL; tmp=tmp->next)
  {
     if (tmp->entityID==eID && tmp->itemID==iID)
        tmp->visible = FALSE;
  }
}

extern void set_seqvisible (Uint2 eID, Uint2 iID, SelEdStructPtr sesp)
{
  SelEdStructPtr tmp;
  for (tmp=sesp;tmp!=NULL; tmp=tmp->next)
  {
     if (tmp->entityID==eID && tmp->itemID==iID)
        tmp->visible = TRUE;
  }
}
extern Boolean is_seqvisible (Uint2 eID, Uint2 iID, SelEdStructPtr sesp)
{
  SelEdStructPtr tmp;
  for (tmp=sesp;tmp!=NULL; tmp=tmp->next)
  {
     if (tmp->entityID==eID && tmp->itemID==iID)
        return (Boolean)tmp->visible;
  }
  return FALSE;
}

/******************************************************************/
extern void setposition_toses (SelEdStructPtr ssp, Int4 from, Int4 to)
{
  SeqLocPtr        slp;
  SeqIntPtr        sitcaret;
  if (ssp == NULL) return;
  if (ssp->region == NULL) return;
  slp = (SeqLocPtr) ssp->region;
  sitcaret = (SeqIntPtr) slp->data.ptrvalue;
  sitcaret->from = from;
  sitcaret->to = to;
}


extern SelEdStructPtr ss_to_ses (SelStructPtr ssp)
{
  SeqLocPtr      slp, slpses;
  SelEdStructPtr sesp;
  SeqIdPtr       sip;
  BioseqPtr      bsp;

  if (ssp == NULL) return NULL;
  sesp = (SelEdStructPtr) MemNew (sizeof (SelEdStruct));
  sesp->entityID = ssp->entityID;
  sesp->itemID = ssp->itemID;
  sesp->itemtype = ssp->itemtype;
  if ( ssp->regiontype == OM_REGION_SEQLOC) {
         sesp->regiontype = OM_REGION_SEQLOC;
         if ( ssp->region == NULL ) {
               ErrPostEx(SEV_WARNING, 0, 0, "Fail in SelStructCpy [1]");
               return NULL;
         }
         slp = (SeqLocPtr) ssp->region;
         sip = NULL;
         bsp = BioseqLockById (SeqLocId (slp));
         if (bsp != NULL) {
            for (sip=bsp->id; sip!= NULL; sip = sip->next) {
               if (sip->choice == SEQID_GI)
                  break; 
            }
            BioseqUnlock (bsp); 
         } 
         if (sip == NULL) {
            sip = SeqLocId (slp);
            sip = SeqIdFindBest (sip, 0);
         }
         slpses = SeqLocIntNew (SeqLocStart(slp), SeqLocStop(slp), SeqLocStrand(slp), sip);
         if ( slpses == NULL ) {
               ErrPostEx(SEV_WARNING, 0, 0, "Fail in SelStructCpy [2]");
               return NULL;
         }
         sesp->region = (Pointer) slpses;
  } 
  else {
         sesp->regiontype = 0;
         sesp->region = NULL;
  }
  sesp->label[0] ='\0';
  sesp->data = NULL;
  sesp->next = NULL;
  sesp->prev = NULL;
  return sesp;
}

extern SelStructPtr ses_to_ss (SelEdStructPtr sesp)
{
  SeqLocPtr    slp, slpses;
  SelStructPtr ssp;

  if (sesp == NULL) return NULL;
  ssp = (SelStructPtr) MemNew (sizeof (SelStruct));
  ssp->entityID = sesp->entityID;
  ssp->itemID = sesp->itemID;
  ssp->itemtype = sesp->itemtype;
  if ( sesp->regiontype == OM_REGION_SEQLOC ) {
         ssp->regiontype = OM_REGION_SEQLOC;
         if ( sesp->region == NULL ) {
               ErrPostEx(SEV_WARNING, 0, 0, "Fail in SelStructCpy [1]");
               return NULL;
         }
         slpses = (SeqLocPtr) sesp->region;
         slp = SeqLocIntNew (SeqLocStart(slpses), SeqLocStop(slpses), 
                              SeqLocStrand(slpses), SeqLocId(slpses));
         if ( slp == NULL ) {
               ErrPostEx(SEV_WARNING, 0, 0, "Fail in SelStructCpy [2]");
               return NULL;
         }
         ssp->region = (Pointer) slp;
  } 
  else {
         ssp->regiontype = 0;
         ssp->region = NULL;
  }
  ssp->next = NULL;
  ssp->prev = NULL;
  return ssp;
}

/******************************************************************/
extern Boolean is_samess_ses (SelStructPtr ssp1, SelEdStructPtr ssp2)
{
  if (ssp1 == NULL || ssp2 == NULL) return FALSE;
  if (ssp1->entityID != ssp2->entityID) return FALSE;
  if (ssp1->itemID != ssp2->itemID) return FALSE;
  if (ssp1->itemtype != ssp2->itemtype) return FALSE;
  return TRUE;
}

/******************************************************************/
extern Boolean is_sameses (SelEdStructPtr ssp1, SelEdStructPtr ssp2)
{
  if (ssp1 == NULL || ssp2 == NULL) return FALSE;
  if (ssp1->entityID != ssp2->entityID) return FALSE;
  if (ssp1->itemID != ssp2->itemID) return FALSE;
  if (ssp1->itemtype != ssp2->itemtype) return FALSE;
  return TRUE;
}

/*********************************************************
***
***  ObjMgr procedures
***
**********************************************************/
extern Boolean AfterAlsoSelect (void)
{
  SelStructPtr sspa = NULL,
               sspb = NULL;
  SeqLocPtr    slpa,
               slpb;
  SeqIntPtr    sint;
  Uint2        eIDa, iIDa, ita,
               eIDb, iIDb, itb;
  Boolean      check = TRUE,
               loopin = TRUE,
               changed = FALSE;
  
  while (check) 
  {
     check = FALSE;
     loopin = TRUE;
     sspa = ObjMgrGetSelected ();
     while (sspa != NULL && loopin) 
     {
        if ( checkssp_for_editor (sspa)) 
        {
           for (sspb = sspa->next; sspb != NULL; sspb = sspb->next) 
           {
              eIDa = sspa->entityID;
              iIDa = sspa->itemID;
              ita = sspa->itemtype;
              eIDb = sspb->entityID;
              iIDb = sspb->itemID;
              itb = sspb->itemtype;
              if ( checkssp_for_editor (sspb) && is_sameId (eIDa, iIDa, ita, 255, eIDb, iIDb, itb, 255) )
              {
                 slpa = (SeqLocPtr)sspa->region;
                 slpb = (SeqLocPtr)sspb->region;
                 if (SeqLocCompare (slpa, slpb) == SLC_A_IN_B) {
                    ObjMgrDeSelect(eIDa, iIDa, ita, sspa->regiontype, slpa);
                    check = TRUE;
                    loopin = FALSE;
                    changed = TRUE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_B_IN_A) {
                    ObjMgrDeSelect(eIDb, iIDb, itb, sspb->regiontype, slpb);
                    check = TRUE;
                    loopin = FALSE;
                    changed = TRUE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_A_EQ_B) {
                    ObjMgrDeSelect(eIDb, iIDb, itb, sspb->regiontype, slpb);
                    check = TRUE;
                    loopin = FALSE;
                    changed = TRUE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_A_OVERLAP_B) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    if (SeqLocStart(slpa) < SeqLocStart(slpb)) {
                       sint->to = SeqLocStop(slpb);
                    }
                    else {
                       sint = (SeqIntPtr) slpa->data.ptrvalue;
                       sint->from = SeqLocStart(slpb);
                    }
                    ObjMgrDeSelect(eIDb, iIDb, itb, sspb->regiontype, slpb);
                    check = TRUE;
                    loopin = FALSE;
                    changed = TRUE;
                    break;
                 }
                 else if (SeqLocStop(slpa) == SeqLocStart(slpb)-1) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    sint->to = SeqLocStop(slpb);
                    ObjMgrDeSelect(eIDb, iIDb, itb, sspb->regiontype, slpb);
                    check = TRUE;
                    loopin = FALSE;
                    changed = TRUE;
                    break;
                 }
                 else if (SeqLocStart(slpa) == SeqLocStop(slpb)+1) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    sint->from = SeqLocStart(slpb);
                    ObjMgrDeSelect(eIDb, iIDb, itb, sspb->regiontype, slpb);
                    check = TRUE;
                    loopin = FALSE;
                    changed = TRUE;
                    break;
                 }
              }
           }
        }
        if (loopin && sspa != NULL)
           sspa = sspa->next;
     }
  }
  return changed;
}

extern void ObjMgrSelectPrint (void)
{
  SelStructPtr ssp = NULL;
  FILE *fout;
  Char    strLog[128];

  fout = FileOpen("LogFile", "a");
  if (fout==NULL) 
     return;

  fprintf (fout, "ObjMgrSelectPrint\n");
  ssp = ObjMgrGetSelected();  
  for (; ssp != NULL; ssp = ssp->next) 
  {
    if (ssp->regiontype == OM_REGION_SEQLOC) {
     SeqIdWrite(SeqLocId((SeqLocPtr)ssp->region),strLog,PRINTID_FASTA_LONG,120);
     if ( ssp->region != NULL ) {
        fprintf (fout, "selstruc  %d %d %d  %s %ld %ld \n", (int)ssp->entityID, (int)ssp->itemID, (int)ssp->itemtype, strLog, (long)SeqLocStart((SeqLocPtr)ssp->region), (long)SeqLocStop((SeqLocPtr)ssp->region));
     }
     else 
        fprintf (fout, "selstruc %d %d %d region=NULL\n", (int)ssp->entityID, (int)ssp->itemID, (int)ssp->itemtype);
    }
    else 
     fprintf (fout, "selstruc %d %d %d regiontype=0\n", (int)ssp->entityID, (int)ssp->itemID, (int)ssp->itemtype);
  }
  fprintf (fout, "\n");
  FileClose(fout);
}

/******************************************************************/
extern void SelectType (EditAlignDataPtr adp, Uint2 feattype, Int4 slpto)
{
  ValNodePtr       vnp;
  SeqLocPtr        slp;
  SelStructPtr     ssp;
  Uint2            ei, ii, it;
  AlignNodePtr     anp;
  Boolean          first = TRUE;
 
  if (adp == NULL) return;
  if (adp->anp_list == NULL) return;
  for (vnp =adp->anp_list; vnp !=NULL; vnp =vnp->next)
  {
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     if (anp != NULL)
     {
           ei = anp->seq_entityID;
           ii = anp->bsp_itemID;
           it = feattype;
           ssp = is_selectedbyID (ei, ii, it);
           if (ssp == NULL)
           {
              slp = SeqLocIntNew (0, slpto, Seq_strand_plus, anp->sip);
              if (first) {
                 ObjMgrSelect (ei, ii, it, OM_REGION_SEQLOC, slp);
                 first = FALSE;
              }
              else
                 ObjMgrAlsoSelect (ei, ii, it, OM_REGION_SEQLOC, slp);
           }
     }
  }
  return;
}

/******************************************************************/
extern Int2 GetNumberObjMgrSelect (void)
{
  SelStructPtr ssp = NULL;
  Int2         nselect = 0;

  ssp = ObjMgrGetSelected ();
  for (; ssp != NULL; ssp = ssp->next) 
     if (checkssp_for_editor(ssp)) nselect++;
  return nselect;
}

/******************************************************************/
extern Int2 checkOMss_for_itemtype (Uint2 itemtype)
{
  SelStructPtr ssp = NULL;
  Int2         nselect = 0;

  ssp = ObjMgrGetSelected ();
  for (; ssp != NULL; ssp = ssp->next) 
     if (ssp->itemtype == itemtype && checkssp_for_editor (ssp)) 
        nselect++;
  return nselect;
}

/******************************************************************/
extern SelStructPtr getOMselect_for_itemtype (Uint2 itemtype)
{
  SelStructPtr ssp = NULL;
  Int2         nselect = 0;

  ssp = ObjMgrGetSelected ();
  for (; ssp != NULL; ssp = ssp->next) 
     if (ssp->itemtype == itemtype && checkssp_for_editor (ssp)) 
        break;
  return ssp;
}

/******************************************************************/
extern SelStructPtr is_selectedbyID (Uint2 entityID, Uint2 itemID, Uint2 itemtype)
{
  SelStructPtr ssp = NULL;

  ssp = ObjMgrGetSelected();  
  for (; ssp != NULL; ssp = ssp->next) 
  {
     if ( is_sameId (ssp->entityID, ssp->itemID, ssp->itemtype, 255, entityID, itemID, itemtype, 255) ) 
     {
        break;
     }
  }
  return ssp;
}

extern SelEdStructPtr getCDSselect (ValNodePtr seqfeathead, ValNodePtr feathead)
{
  SelStructPtr   ssp = NULL;
  SelEdStructPtr feat = NULL;

  ssp = ObjMgrGetSelected ();
  for (; ssp != NULL; ssp = ssp->next) if (checkssp_for_editor (ssp)) {
     if ( ssp->itemtype == OBJ_SEQFEAT ) 
     {
        feat = get_feat_fromid (seqfeathead, FEATDEF_CDS, ssp->entityID, ssp->itemID, -1, NULL);
        if (feat != NULL) break;
     }
     else if ( ssp->itemtype == OBJ_VIRT ) 
     {
        feat = get_feat_fromid (feathead, SEQFEAT_CDREGION, ssp->entityID, ssp->itemID, -1, NULL);
        if (feat != NULL) break;
     }
  }
  return feat;
}

extern Int2 checkCDSselect_forprotein (ValNodePtr seqfeathead, ValNodePtr feathead, Boolean with_prot)
{
  SelStructPtr   ssp = NULL;
  SelEdStructPtr feat = NULL;
  Int2           nselect = 0;

  ssp = ObjMgrGetSelected ();
  for (; ssp != NULL; ssp = ssp->next) if (checkssp_for_editor (ssp)) {
     if ( ssp->itemtype == OBJ_SEQFEAT ) 
     {
        feat = get_feat_fromid (seqfeathead, FEATDEF_CDS, ssp->entityID, ssp->itemID, -1, NULL);
        if (feat != NULL)
           if (with_prot && feat->data != NULL) nselect++;
           else if (!with_prot && feat->data == NULL) nselect++;
     }
     else if ( ssp->itemtype == OBJ_VIRT ) 
     {
        feat = get_feat_fromid (feathead, SEQFEAT_CDREGION, ssp->entityID, ssp->itemID, -1, NULL);
        if (feat != NULL)
           if (with_prot && feat->data != NULL) 
              nselect++;
           else if (!with_prot && feat->data == NULL) 
              nselect++;
     }
  }
  return nselect;
}


/******************************************************************
***   checkssp_for_editor
******************************************************************/
extern Boolean checkssp_for_editor (SelStructPtr ssp)
{
  SeqLocPtr slp;

  if (ssp == NULL) return FALSE;
  if (ssp->regiontype == 0) return FALSE;
  if (ssp->region == NULL) {
         ssp->regiontype = 0;
         return FALSE;
  }
  slp = (SeqLocPtr)ssp->region;
  if (SeqLocStart(slp)<0) {
  }
  if (SeqLocStop(slp)==-2) {
  }
  return TRUE;
}

/***********************************************************
***
***
************************************************************/
extern SeqLocPtr checkseqlocfeature_for_editor (Uint2 entityID, Uint2 itemID, ValNodePtr headfeat)
{
  SelEdStructPtr feat;
  SeqLocPtr      slp = NULL,
                 tmp;
  SeqIdPtr       sip;
  Int4           start, stop;
  Uint1          strand;

  if (headfeat == NULL) 
     return NULL;
  feat = get_feat_fromid (headfeat, 255, entityID, itemID, -1, NULL);
  if (feat != NULL) 
  {
     tmp = (SeqLocPtr) feat->region;
     sip = SeqLocId (tmp);
     strand = SeqLocStrand (tmp);
     start = SeqLocStart (tmp);
     if (feat->next != NULL)
        while (feat->next != NULL) feat = feat->next;
     stop = SeqLocStop ((SeqLocPtr) feat->region);
     slp = SeqLocIntNew (start, stop, strand, sip);
  }
  return slp;
}

extern void checkselectsequinfeature_for_editor (ValNodePtr headfeat)
{
  SelStructPtr   ssp;
  SelEdStructPtr feat;
  SeqLocPtr      slp, tmp;
  SeqIdPtr       sip;
  Int4           start, stop;
  Int2           k = 0;
  Uint1          strand;

  if (headfeat == NULL) return;
  for (ssp = ObjMgrGetSelected (); ssp != NULL; ssp = ssp->next) 
          if (ssp->itemtype == OBJ_SEQFEAT) k++;
  if (k > 0) {
     ssp = ObjMgrGetSelected ();
     for (; ssp != NULL; ssp = ssp->next) {
        if (ssp->itemtype == OBJ_SEQFEAT) 
        {
           if (ssp->regiontype != OM_REGION_SEQLOC || ssp->region == NULL)
           {
              feat = get_feat_fromid (headfeat, 255, ssp->entityID, ssp->itemID, -1, NULL);
              if (feat != NULL) {
                 tmp = (SeqLocPtr) feat->region;
                 sip = SeqLocId (tmp);
                 strand = SeqLocStrand (tmp);
                 start = SeqLocStart (tmp);
                 if (feat->next != NULL)
                    while (feat->next != NULL) feat = feat->next;
                 stop = SeqLocStop ((SeqLocPtr) feat->region);
                 slp = SeqLocIntNew (start, stop, strand, sip);
                 ssp->regiontype = OM_REGION_SEQLOC;
                 ssp->region = slp;
              }
           }
        } 
     }
  }
}

extern Int4 getminpos_fromOMselect (Uint2 itemsubtype)
{
  SelStructPtr   ssp;
  Int4           minpos = INT4_MAX;
  Int4           from;

  ssp = ObjMgrGetSelected ();
  for (; ssp != NULL; ssp = ssp->next) 
  {
     if (ssp->itemtype == itemsubtype) 
     {
        if (ssp->regiontype == OM_REGION_SEQLOC && ssp->region != NULL)
        {
            from = SeqLocStart ((SeqLocPtr)ssp->region);
            if (from < minpos) minpos = from;
        } 
     }
  }
  if ( minpos < INT4_MAX ) return minpos;
  return -1;
}

/***********************************************
***  locate_in_seqalign   
***    in : pos in Align coordinates
***    out: seen TRUE if pos in salp
************************************************/
extern Boolean locate_in_seqalign (Int4 pos, Int2 dim, Int2 dspnumseg, BoolPtr *dspstart, Int4Ptr *dsplens, Int2 *numseg_before, Int2 *subdsplens, Int4 *sumdsplens_before)
{
  BoolPtr     start = *dspstart;
  Int4Ptr     lens  = *dsplens;
  Int4        sumlens= 0;
  Int4        sumlensseq= 0;
  Int2        numseg = 0;
  Int2        sublens;
  Boolean     seen = FALSE;

  if ( dspnumseg == 0 || start == NULL || lens == NULL ) {
         ErrPostEx(SEV_WARNING, 0, 0, "fail in locate_in_seqalign [1]\n");
         return FALSE;
  }
  while ( !seen && numseg < dspnumseg ) {
         numseg++;
         if ( pos  >= sumlens && pos < sumlens + *lens ) {
                sublens = abs (pos - sumlens);
                seen = TRUE;
         }
         else {
                if ((Boolean)(*start)) sumlensseq += *lens;
                if ( numseg == dspnumseg ) break;
                start += dim; 
                sumlens += *lens;
                lens++;
         }
  }
  if ( seen )
  {
        *dspstart  = start;
        *dsplens   = lens;
        *numseg_before = numseg;
        *subdsplens= sublens;
        *sumdsplens_before= sumlensseq;
  }
  return seen;
}

/************************************
*** SeqCoordToAlignCoord
**
************************************/
extern Int4 SeqCoordToAlignCoord (Int4 position, SeqIdPtr sip, SeqAlignPtr salp, Int2 intersalpwidth, Int2 is_end)
{
  CompSegPtr  dsp;
  BoolPtr     dspstart;
  Int4Ptr     dsplens;
  Int4        from;
  Int4        sumlens = 0;
  Int4        seqlens = 0;
  Int4        lensplus;
  Int4        start, stop;
  Int2        numseg;
  Int2        inter_salp = 0;
  Int2        index;
  Uint1       dspstrand = Seq_strand_unknown;
  Boolean     seen = FALSE;

  if (is_end == NO_RESIDUE)
     return position;
  if (position < 0)
     return position;

  dsp = (CompSegPtr) salp->segs;
  if (dsp == NULL) {
     return GAP_RESIDUE;
  }
  index = position_inIdlist (sip, dsp->ids);
  if (index == 0) {
     return GAP_RESIDUE;
  }
  index -= 1;
  from = *(dsp->from + index);
  if (is_end == GAP_RESIDUE)
     position = from;
  dspstart = dsp->starts + index;
  dsplens = dsp->lens;
  if (dspstart == NULL || dsplens == NULL ) {
     return GAP_RESIDUE;
  }
  if (dsp->strands!=NULL)
     dspstrand = *(dsp->strands + index);
  
  numseg = 1;
  while ( !seen && numseg <= dsp->numseg ) 
  {
     if (dspstrand ==Seq_strand_minus) {
        start= from - seqlens - *dsplens;
        stop = from - seqlens;
     } else {
        start= from + seqlens;
        stop = from + seqlens + *dsplens;
     }
     if (*dspstart && position >= start && position < stop) 
     {
        if (dspstrand ==Seq_strand_minus)
           lensplus = abs (from + seqlens - position);
        else
           lensplus = abs (position - from - seqlens);
        seen = TRUE;
     }
     else if (*dspstart && position <= stop) {
/**
        if (is_end == APPEND_RESIDUE ) 
**/
        {
           if (dspstrand ==Seq_strand_minus)
              lensplus = abs (from + seqlens - position);
           else
              lensplus = abs (position - from - seqlens);
           seen = TRUE;
        }
     }
     else if ( numseg == dsp->numseg ) 
     {
        if ( salp->next == NULL ) break; 
        else 
        { 
           salp = salp->next;
           dsp = (CompSegPtr) salp->segs;
           from = *(dsp->from + index);
           dspstart = dsp->starts + index;
           dsplens = dsp->lens;
           inter_salp++;
           numseg = 1;
        }
     }
     else if (numseg < dsp->numseg) 
     {
           sumlens += *dsplens;
           if (*dspstart) 
              seqlens += *dsplens;
           dspstart += dsp->dim; 
           dsplens++;
     }
     numseg++;
  }
  if ( !seen ) {
     if (!(*dspstart))    /***** if after sequence 2 mais seq1 last segment***/
         return seqlens; 
     return GAP_RESIDUE;
  }
  if (position == APPEND_RESIDUE)
     return position;
  position = sumlens + lensplus + intersalpwidth*inter_salp;
  return position;
}

/************************************
*** AlignCoordToSeqCoord
************************************/
extern Int4 AlignCoordToSeqCoord (Int4 position, SeqIdPtr sip, SeqAlignPtr salp,ValNodePtr sqloc_list, Int2 intersalpwidth)
{
  CompSegPtr  dsp;
  BoolPtr     dspstart;
  Int4Ptr     dsplens;
  Int4        from;
  Int4        sumlens = 0;
  Int4        sumstart = 0;
  Int4        seqlens = 0;
  Int2        numseg = 0;
  Int2        inter_salp = 0;
  Int2        index;
  Uint1       dspstrand = Seq_strand_unknown;
  Boolean     seen = FALSE;

  if (position == APPEND_RESIDUE)
     return position;
  dsp = (CompSegPtr) salp->segs;
  if (dsp == NULL) 
     return (Int4)GAP_RESIDUE;
  index = position_inIdlist (sip, dsp->ids);
  if (index == 0) 
     return (Int4)GAP_RESIDUE;
  index -= 1;
  from = *(dsp->from + index);
/*
  if (position <= from) {
     return from;
  }
*/
  dspstart = dsp->starts + index;
  dsplens = dsp->lens;
  if (dspstart == NULL || dsplens == NULL ) {
     return (Int4)GAP_RESIDUE;
  }
  if (!(*dspstart) && (position < *dsplens)) {
     return (Int4)GAP_RESIDUE;
  }
  if (dsp->strands!=NULL)
     dspstrand = *(dsp->strands + index);
  numseg = 0;
  sumlens = 0;
  for (; dsplens != NULL && numseg < dsp->numseg; dsplens++, numseg++) 
     sumlens += *dsplens;
  dsplens = dsp->lens; 
  if (position >= sumlens) {     
     sumlens = 0;
     numseg = 0;
     for (; dsplens != NULL && numseg < dsp->numseg; dsplens++, numseg++) {
        if (*dspstart)  sumlens += *dsplens;
        dspstart += dsp->dim; 
     }
    seen = TRUE;
     if (dspstrand == Seq_strand_minus) {
        position = from - sumlens;
     } else {
        position = from + sumlens;
     }
  }
  else {
        sumlens = 0;
        numseg = 0;
        while ( !seen && numseg < dsp->numseg ) {
            numseg++;
            if (position >=sumlens && position <sumlens +*dsplens ) {
               if (*dspstart) 
                  seqlens += abs (position - sumlens);
               seen = TRUE;
            }
            else if ( numseg == dsp->numseg ) 
            {
               if ( salp->next == NULL ) break; 
               else 
               { 
                  sumstart += sumlens + *dsplens;
                  salp = salp->next;
                  dsp = (CompSegPtr) salp->segs;
                  from = *(dsp->from + index);
                  dspstart = dsp->starts + index;
                  dsplens = dsp->lens;
                  inter_salp++;
                  numseg = 0;
               }
            }
            else 
            {
               if ( *dspstart ) 
                  seqlens += *dsplens;
               sumlens += *dsplens;
               dspstart += dsp->dim; 
               dsplens++;
            }
        }
        if (seen) {
           if (dspstrand == Seq_strand_minus) {
              position = from - seqlens - sumstart;
           } else {
              position = from + seqlens - sumstart;
           }
     }
  }
  if ( !seen ) { 
     return (Int4)GAP_RESIDUE;
  }
  index =chkloc (sip, position, sqloc_list, &from);
  if (index == GAP_RESIDUE || index == APPEND_RESIDUE) {
     return (Int4)index;
  }
  return position;
}

extern Int4 AlignCoordToSeqCoord2 (Int4 position, SeqIdPtr sip, SeqAlignPtr salp,ValNodePtr sqloc_list, Int2 intersalpwidth)
{
  CompSegPtr  dsp;
  BoolPtr     dspstart;
  Int4Ptr     dsplens;
  Int4        from;
  Int4        sumlens = 0;
  Int4        sumstart = 0;
  Int4        seqlens = 0;
  Int2        numseg = 0;
  Int2        inter_salp = 0;
  Int2        index;
  Uint1       dspstrand = Seq_strand_unknown;
  Boolean     seen = FALSE;

  if (position == APPEND_RESIDUE)
     return position;
  dsp = (CompSegPtr) salp->segs;
  if (dsp == NULL) 
     return (Int4)GAP_RESIDUE;
  index = position_inIdlist (sip, dsp->ids);
  if (index == 0) 
     return (Int4)GAP_RESIDUE;
  index -= 1;
  from = *(dsp->from + index);
/*
  if (position <= from) {
     return from;
  }
*/
  dspstart = dsp->starts + index;
  dsplens = dsp->lens;
  if (dspstart == NULL || dsplens == NULL ) {
     return (Int4)GAP_RESIDUE;
  }
  if (!(*dspstart) && (position < *dsplens)) {
     return (Int4)GAP_RESIDUE;
  }
  if (dsp->strands!=NULL)
     dspstrand = *(dsp->strands + index);
  numseg = 0;
  sumlens = 0;
  for (; dsplens != NULL && numseg < dsp->numseg; dsplens++, numseg++) 
     sumlens += *dsplens;
  dsplens = dsp->lens; 
  if (position >= sumlens) {     
     position = (Int4)APPEND_RESIDUE;
/*****
     sumlens = 0;
     numseg = 0;
     for (; dsplens != NULL && numseg < dsp->numseg; dsplens++, numseg++) {
        if (*dspstart)  sumlens += *dsplens;
        dspstart += dsp->dim; 
     }
    seen = TRUE;
     if (dspstrand == Seq_strand_minus) {
        position = from - sumlens;
     } else {
        position = from + sumlens;
     }
****/
  }
  else {
        sumlens = 0;
        numseg = 0;
        while ( !seen && numseg < dsp->numseg ) {
            numseg++;
            if (position >=sumlens && position <sumlens +*dsplens ) {
               if (*dspstart) 
                  seqlens += abs (position - sumlens);
               seen = TRUE;
            }
            else if ( numseg == dsp->numseg ) 
            {
               if ( salp->next == NULL ) break; 
               else 
               { 
                  sumstart += sumlens + *dsplens;
                  salp = salp->next;
                  dsp = (CompSegPtr) salp->segs;
                  from = *(dsp->from + index);
                  dspstart = dsp->starts + index;
                  dsplens = dsp->lens;
                  inter_salp++;
                  numseg = 0;
               }
            }
            else 
            {
               if ( *dspstart ) 
                  seqlens += *dsplens;
               sumlens += *dsplens;
               dspstart += dsp->dim; 
               dsplens++;
            }
        }
        if (seen) {
           if (dspstrand == Seq_strand_minus) {
              position = from - seqlens - sumstart;
           } else {
              if (!(*dspstart)) {
                 position =chkloc (sip, position, sqloc_list, &from);
                 if (position==0) 
                    position = (Int4)GAP_RESIDUE;
              }
              else
                 position = from + seqlens - sumstart;
           }
     }
  }
  if ( !seen ) { 
     position =chkloc (sip, position, sqloc_list, &from);
  }
  return position;
}


extern Boolean GetAlignCoordFromSeqLoc (SeqLocPtr slp, SeqAlignPtr salp, Int4 *start, Int4 *stop)
{
  SeqIdPtr         sip;
  Int4             from, to;
 
  sip = SeqLocId (slp);
  from = SeqLocStart (slp);
  from = SeqCoordToAlignCoord (from, sip, salp, 0, 0);
  to = SeqLocStop (slp);
  to = SeqCoordToAlignCoord (to, sip, salp, 0, 0);
  if (from < 0 || to < 0)
     return FALSE;
  *start = from;
  *stop = to;
  return TRUE;
}
 
/************************************************************************
***
***   SeqPosToLineColumn
***
************************************************************************/
extern Boolean SeqPosToLineColumn (Uint2 itemID, Uint2 entityID, Uint2 itemtype, Int4 pos, Int2 *line, Int2 *column, Int4 hoffset, EditAlignDataPtr adp)
{
  Boolean seen = FALSE;
  Int4    from, to;
  Int2    col = 0;
  Int4    length;
  Int2    j;

  *line = -1;
  *column = -1;
  if ( pos < 0 ) 
         return FALSE;
  hoffset -= adp->bufferstart;
  length = adp->length - adp->bufferstart;
  for ( j=0; !seen && j < MAXLineWindow; j++ )
  {
         if (itemID == adp->item_id[j] && entityID == adp->seqEntity_id[j] && itemtype == adp->itemtype[j] ) 
         { 
                from = hoffset + adp->alignline [j] * adp->visibleWidth;
                if (from > length) break;
                to = hoffset + (adp->alignline [j]+1) * adp->visibleWidth -1;
                if (to > length) to = length;
/**
                WriteLog ("SeqPosToLineColumn %d %ld %ld %ld %d %d %d %d\n", j,  
                       pos, from, to, adp->item_id[j], adp->seqEntity_id[j], 
                       adp->colonne[from], adp->colonne[to]);
**/
                if( pos >= adp->colonne[from]  && pos <= adp->colonne[to]) 
                {
                       col = (Int2)(pos - adp->colonne[from]);
                       if (adp->columnpcell > 0)
                          col += (Int2)(col) / (Int2)(adp->columnpcell);
                       seen = TRUE;
                       break;
                }
         }
  }
  if ( !seen || j == MAXLineWindow ) 
         return FALSE;
  *line = j;
  *column = col;
  return TRUE;
}

/************************************************************************
***
***   SeqPosInLineColumn
***
************************************************************************/
extern Boolean SeqPosInLineColumn (Int4 pos, Int2 alignline, Int2 *column, Int4 hoffset, EditAlignDataPtr adp)
{
  Int4    from, to;
  Int4    length;
  Int4    col = 0;

  if ( pos < 0 )  return FALSE;
  hoffset -= adp->bufferstart;
  length = adp->length - adp->bufferstart;
  from = hoffset + alignline * adp->visibleWidth;
  to = hoffset + (alignline +1) * adp->visibleWidth - 1; 
  if (from > length) return FALSE;
  if (to > length) to = length;

  if ( pos < adp->colonne[from] ) 
  {
         *column = -1;
         return FALSE;
  }
  if ( pos > adp->colonne[to] ) 
  {
         *column = -2;
         return FALSE;
  }
  col = pos - adp->colonne[from];
  if (adp->columnpcell > 0)
         col += (Int4) col / (Int4) adp->columnpcell;
  *column = (Int2) col;
  return TRUE;
}
 
/***********************************************
***  LocateInSeqAlign
***    in : pos in Align coordinates
***    out: seen TRUE if pos in salp
************************************************/
extern Boolean LocateInSeqAlign (Int4 pos, Int2 dim, Int2 dspnumseg, BoolPtr *dspstart, Int4Ptr *dsplens, Int2 *numseg_before, Int4 *subdsplens, Int4 *sumdsplens_before)
{
  BoolPtr     start = *dspstart;
  Int4Ptr     lens  = *dsplens;
  Int4        sumlens= 0;
  Int4        sumlensseq= 0;
  Int2        numseg = 0;
  Int4        sublens;
  Boolean     seen = FALSE;
 
  if ( dspnumseg == 0 || start == NULL || lens == NULL ) {
         ErrPostEx(SEV_WARNING, 0, 0, "fail in LocateInSeqAlign [1]\n");
         return FALSE;
  }
  while ( !seen && numseg < dspnumseg ) {
         numseg++;
         if ( pos  >= sumlens && pos < sumlens + *lens ) {
                sublens = abs (pos - sumlens);
                seen = TRUE;
         }
         else {
                if ((Boolean)(*start)) sumlensseq += *lens;
                if ( numseg == dspnumseg ) break;
                start += dim;
                sumlens += *lens;
                lens++;
         }
  }
  if ( seen )
  {
        *dspstart  = start;
        *dsplens   = lens;
        *numseg_before = numseg;
        *subdsplens= sublens;
        *sumdsplens_before= sumlensseq;
  }
  return seen;
}


/***********************************************
***  LocateInSeqAlign
***    in : pos in Align coordinates
***    out: seen TRUE if pos in salp
************************************************/
extern Boolean LocateInSeqAlignDenSeg (Int4 pos, Int2 dim, Int2 dspnumseg, Int4Ptr *dspstart, Int4Ptr *dsplens, Int2 *numseg_before, Int4 *subdsplens)
{
  Int4Ptr     start = *dspstart;
  Int4Ptr     lens  = *dsplens;
  Int4        sumlens= 0;
  Int2        numseg = 0;
  Int4        sublens= 0;
  Boolean     seen = FALSE;
 
  if ( dspnumseg == 0 || start == NULL || lens == NULL ) {
         ErrPostEx(SEV_WARNING, 0, 0, "fail in LocateInSeqAlignDenSeg [1]");
         return FALSE;
  }
  if (pos == -1 || *start == -1) 
  {
     numseg = 1;      /*!!!!!!!!!!!!!!!!!!!!!*/
     seen = TRUE; 
  }
  else {
     sumlens = *start;
     while ( !seen && numseg < dspnumseg ) {
         numseg++;
         if ( pos  >= sumlens && pos < sumlens + *lens ) {
                sublens = abs (pos - sumlens);
                seen = TRUE;
         }
         else {
                if ( numseg == dspnumseg ) break;
                start += dim;
                sumlens += *lens;
                lens++;
         }
     }
  }
  if ( seen )
  {
        *dspstart  = start;
        *dsplens   = lens;
        *numseg_before = numseg;
        *subdsplens= sublens;
  }
  return seen;
}
 
/*******************************************************
***
***   FIND functions
***
*******************************************************/
extern SelStructPtr SetupMatchPatList (ValNodePtr match_pat, ValNodePtr anp_list)
{
  SelStructPtr    head = NULL;
  SelStructPtr    ssp,
                  tmp;
  SeqLocPtr       slp;
  AlignNodePtr    anp;
  ValNodePtr      vnp, vnp2;

  if (match_pat == NULL || anp_list == NULL)
     return NULL;
  for (vnp2 = match_pat; vnp2!=NULL; vnp2=vnp2->next)
  {
     slp = (SeqLocPtr) vnp2->data.ptrvalue;
     if (slp != NULL) {
        for (vnp = anp_list; vnp != NULL; vnp = vnp->next) {
           anp = (AlignNodePtr) vnp->data.ptrvalue;
           if (SeqIdForSameBioseq(anp->sip, SeqLocId(slp)))
              break;
        }
        if (vnp != NULL) 
        {
           ssp = SelStructNew (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ, SeqLocStart(slp), SeqLocStop(slp), SeqLocId(slp), Seq_strand_plus, FALSE);
           if (head == NULL)
              head = ssp;
           else {
              for (tmp=head; tmp->next!=NULL; tmp=tmp->next)
                 continue;
              tmp->next = ssp;
              ssp->prev = tmp;
           }
        }
     }
  }
  ValNodeFreeType (&match_pat, TypeSeqLoc);
  return head;
}

extern SelStructPtr ShowNextPattern (SelStructPtr match_pat, SelStructPtr cur_pat, Int4 start)
{
  SeqLocPtr     slp;
  SelStructPtr  cur;
  Int4          from;
 
  if (match_pat == NULL) 
     return NULL; 
  if (cur_pat == NULL) {
     cur_pat = match_pat;
  }
  else if (cur_pat->next == NULL) {
     if (cur_pat == match_pat) 
        Beep();
     else cur_pat = match_pat;
  }
  else {
     cur_pat = cur_pat->next;
  }
  if (start > 0 && cur_pat != NULL) {
     from = SeqLocStart((SeqLocPtr)cur_pat->region);
     if (from < start) {
        for (cur = cur_pat; cur != NULL; cur = cur->next) {
           from = SeqLocStart((SeqLocPtr)cur->region);
           if (from > start)
              break;
        }
        if (cur != NULL)
           cur_pat = cur;
     }
  }
  if (cur_pat != NULL) {
     if (cur_pat->regiontype == OM_REGION_SEQLOC && cur_pat->region != NULL) {
        slp = AsnIoMemCopy (cur_pat->region, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ObjMgrSelect(cur_pat->entityID, cur_pat->itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, (Pointer)slp);
     }
  }
  return cur_pat;
}

extern SelStructPtr ShowPrecPattern (SelStructPtr match_pat, SelStructPtr cur_pat, Int4 start)
{
  SelStructPtr  tmp;
  SeqLocPtr     slp;
 
  if (match_pat == NULL) 
     return NULL; 
  if (cur_pat == NULL) {
     cur_pat = match_pat;
  }
  else if (cur_pat->prev == NULL) {
     for (tmp=match_pat; tmp->next!=NULL; tmp=tmp->next) 
        continue;
     cur_pat = tmp;
     if (cur_pat == match_pat) 
        Beep();
  }
  else {
     cur_pat = cur_pat->prev;
  }
  if (cur_pat != NULL) {
     if (cur_pat->regiontype == OM_REGION_SEQLOC && cur_pat->region != NULL) {
        slp = AsnIoMemCopy (cur_pat->region, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ObjMgrSelect(cur_pat->entityID, cur_pat->itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, (Pointer)slp);
     }
  }
  return cur_pat;
}

/**************************************
***  
***  EDITOR Procedure
***  
***************************************/
static Char insertstr[] = "Insertboard";

extern Boolean insertchar (CharPtr str, Int4 pos, SeqIdPtr target, Uint1 mol_type, Boolean spliteditmode)
{
  SeqEntryPtr      sep;
  Int4             strlens;
  SeqId            insertboard;
  ObjectId         oip;
  Boolean          ok = FALSE;

  if ( str == NULL ) 
     return FALSE;
  strlens = StringLen (str);
  if (strlens < 1) 
     return FALSE;
  if (pos < -2 || target == NULL) 
     return FALSE; 
  insertboard.choice = SEQID_LOCAL;
  oip.str = insertstr;
  insertboard.data.ptrvalue = (ObjectIdPtr) &oip; 
  sep = StringToSeqEntry (str, &insertboard, strlens, mol_type);
  if (sep != NULL) {
     ok = BioseqInsert(&insertboard, FIRST_RESIDUE, LAST_RESIDUE, Seq_strand_plus, target, pos, TRUE, TRUE, spliteditmode);
     SeqEntryFree (sep);
  }
  else 
     ErrPostEx(SEV_WARNING, 0, 0,  "fail in insertchar [2]");
  return ok;
}

extern Boolean insertchar_atcaret (CharPtr str, EditAlignDataPtr adp)
{
  SeqLocPtr        slpcaret;
  Int4             pos;
  Int4             strlens;

  if ( str == NULL ) 
     return FALSE;
  strlens = StringLen (str);
  if (strlens < 1) 
     return FALSE;
  if ( adp->caret.regiontype == 0 ) 
     return FALSE;
  slpcaret = (SeqLocPtr) adp->caret.region;
  pos = SeqLocStart (slpcaret); 
  if ( pos == adp->length ) 
     pos = -2;
  return (Boolean) insertchar (str, pos, SeqLocId (slpcaret), adp->mol_type, adp->spliteditmode);
}

extern CharPtr char_to_insert (Char *ch, Int4 lens, Uint1 mol_type)
{
  CharPtr      str, strp;
  Int4         j;

  str = MemNew ((size_t)((lens+2)*sizeof(Char)));
  for (j=0, strp = ch; j<lens && *strp!='\0'; j++, strp++)
     str[j] = *strp;
  str[j] = '\0';
  lens = StringLen (str);

  if ( ISA_na (mol_type)) 
  {
     for (j=0, strp = str; j<lens; j++, strp++) {
        *strp = TO_LOWER(*strp);
        if ( StringChr ("abcdghkmnrstuvwy", *strp) == NULL ) {
           MemFree (str);
           return NULL;
        }
     }
  }
  else if ( ISA_aa (mol_type) )
  {
     for (j=0, strp = str; j<lens; j++, strp++) {
        *strp = TO_UPPER(*strp);
        if ( StringChr ("ABCDEFGHIKLMNPQRSTUVWXYZ-*", *strp) == NULL ) {
           MemFree (str);
           return NULL;
        }
     }
  }
  else {
     for (j=0, strp = str; j<lens; j++, strp++) {
        *strp = TO_UPPER(*strp);
     }
  }
  return str;
}


extern SeqEntryPtr getfirst_sep(SeqEntryPtr sep, Uint1 bsp_mol)
{
  SeqEntryPtr  septmp = NULL;
  SeqEntryPtr  sep2 = NULL;
  BioseqSetPtr bssp;
  BioseqPtr    bsp;
 
  if (sep == NULL) return NULL;
  if (IS_Bioseq_set(sep)) {
         bssp = sep->data.ptrvalue;
         if (bssp == NULL)
                return NULL;
         septmp = (SeqEntryPtr) bssp->seq_set;
         if (septmp == NULL)
                return NULL;
  }
  else if (IS_Bioseq(sep))
     septmp = sep;
 
  while ( septmp != NULL ) {
     if  (IS_Bioseq_set(septmp)) {
        bssp = (BioseqSetPtr) septmp->data.ptrvalue;
        for (sep2 = bssp->seq_set; sep2 != NULL; sep2 = sep2->next) {
           if (IS_Bioseq(sep2)) {
              bsp = (BioseqPtr) sep2->data.ptrvalue;
              if (bsp!=NULL && bsp->id!=NULL) {
                 if (ISA_na(bsp->mol) == ISA_na(bsp_mol)) {
                    return sep2;
                 }
              }   
           }   
        }
     }   
     else if (IS_Bioseq(septmp)) {
        bsp = (BioseqPtr) septmp->data.ptrvalue;
        if (bsp!=NULL && bsp->id!=NULL)
        {
           if (ISA_na(bsp->mol) == ISA_na(bsp_mol)) {
              return septmp;
           }
        }
     }
     septmp = septmp->next;
  }   
  return NULL;
}


/*
static SeqAnnotPtr SeqAnnotForHistSeqAlign (SeqAlignPtr salp)
{
  SeqAnnotPtr      sap = NULL;
  ObjectIdPtr      oip;
  UserFieldPtr     ufp;
  UserObjectPtr    uop;

  sap = SeqAnnotForSeqAlign (salp);
  if (sap != NULL && sap->type == 2) 
  {
          oip = ObjectIdNew ();
          oip->str = StringSave ("Hist Seqalign");
          uop = UserObjectNew ();
          uop->type = oip;

          oip = ObjectIdNew();
          oip->str = StringSave ("Hist Seqalign");
          ufp = UserFieldNew ();
          ufp->choice = 4;
          ufp->data.boolvalue = TRUE;
          ufp->label = oip;

          uop->data = ufp;
          ValNodeAddPointer (&(sap->desc), Annot_descr_user, (Pointer) uop);
  }
  return sap;
}

static void AttachSeqAligntoBioseq (Uint2 entityID, SeqAlignPtr salp)
{
  SeqEntryPtr      sep = NULL;
  SeqAnnotPtr      sap = NULL, tmp;
  BioseqSetPtr     bssp= NULL;
  BioseqPtr        bsp = NULL;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep!=NULL) {
     if (IS_Bioseq(sep))
        bsp = (BioseqPtr) sep->data.ptrvalue;
     else if (IS_Bioseq_set(sep)) { 
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        while (bssp!=NULL && bssp->_class == 7) {
           sep = bssp->seq_set;
           bssp = NULL;
           if (IS_Bioseq(sep))  {
              bsp = (BioseqPtr) sep->data.ptrvalue;  
              break;
           }
           else if (IS_Bioseq_set(sep))
              bssp = (BioseqSetPtr) sep->data.ptrvalue;
        }
     }
     if (bssp!=NULL) {
        if (bssp->annot==NULL) {
           sap = SeqAnnotForHistSeqAlign (salp);
           bssp->annot = sap;
        } else {
           for (tmp=bssp->annot;tmp->next!= NULL; tmp=tmp->next) 
              continue;
           sap = SeqAnnotForHistSeqAlign (salp);
           tmp->next = sap;
        }
     }
     else if (bsp != NULL) { 
        if (bsp->annot==NULL) {
           sap = SeqAnnotForHistSeqAlign (salp);
           bsp->annot = sap;  
        } else { 
           for (tmp=bsp->annot;tmp->next!= NULL; tmp=tmp->next)
              continue;   
           sap = SeqAnnotForHistSeqAlign (salp);
           tmp->next = sap;
        }
     }  
  }
}
*/
/*  
static Uint2 GetEntityIDForBioseqSet (void)
{
  SeqEntryPtr  sep;  
  ValNodePtr   vnp;
  AlignNodePtr anp;
  Uint2        eID = 0;

  sep=GetBestTopParentForItemID (master->entityID, master->itemID, OBJ_BIOSEQ);
  sep = GetTopSeqEntryForEntityID (master->entityID);
  if (sep == NULL) return 0;
  if (IS_Bioseq_set(sep)) {
     eID = SeqMgrGetEntityIDForSeqEntry(sep);
  }
  else {
     for (vnp=anp_list; vnp!=NULL; vnp=vnp->next) {
        anp = (AlignNodePtr) vnp->data.ptrvalue;
        sep = GetBestTopParentForItemID (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ);
        sep = GetTopSeqEntryForEntityID (anp->seq_entityID);
        if (sep!=NULL)
           if (IS_Bioseq_set(sep)) {
              eID = SeqMgrGetEntityIDForSeqEntry(sep);
              break;
           }
     }
  }
  return eID;
}
*/

extern MashPtr MashNew (Boolean is_prot)
{
  MashPtr msp;
  
    if((msp = (MashPtr)MemNew(sizeof(Mash)))==NULL) {
      ErrPostEx(SEV_WARNING, 0, 0, "MemNew returned NULL");
      return NULL;
    }
     msp->is_prot = is_prot;
     msp->band_method = 2;
     msp->multidim = FALSE;
     if (msp->is_prot) 
     {
        msp->reward = 1;
        msp->penalty = -2;
        msp->gap_open = 11;
        msp->gap_extend = 1;
	msp->matrixname="BLOSUM62";
        msp->wordsize = 3;
        msp->gap_x_dropoff = 15;
        msp->gap_x_dropoff_final = 500;   /****25****/
        msp->filter = FILTER_SEG;
     } 
     else {
        msp->reward = 2;
        msp->penalty = -3;
        msp->gap_open = 10;
        msp->gap_extend = 2;
	msp->matrixname=NULL;
        msp->wordsize = 11;
        msp->gap_x_dropoff = 50;
        msp->gap_x_dropoff_final = 50;
        msp->filter = FILTER_DUST;
     }
     msp->translate_prot = FALSE;
     msp->transalp = NULL;
     msp->use_gapped_blast = TRUE;
     msp->lg1_ext = 0; 
     msp->rg1_ext =0; 
     msp->lg2_ext =0; 
     msp->rg2_ext =0; 
     msp->lg1_open =0; 
     msp->lg2_open =0; 
     msp->rg1_open =0; 
     msp->rg2_open =0; 
     msp->blast_threshold = 80;
     msp->choice_blastfilter = 2;      /* 1 */
     msp->splicing = TRUE;              /*FALSE;*/
     msp->map_align = FALSE;
     return msp;
}

static Int4Ptr PNTR LIBCALL BlastStyleMatDestruct(Int4Ptr PNTR matrix){
   matrix=MemFree(matrix);
   return matrix;
}

static GlobalBandStructPtr LIBCALL DestructBandStruct(GlobalBandStructPtr gbsp){
   gbsp->matrix = BlastStyleMatDestruct(gbsp->matrix);
   if(gbsp->seq2!=NULL) gbsp->seq2=MemFree(gbsp->seq2);
   if(gbsp->seq1!=NULL) gbsp->seq1=MemFree(gbsp->seq1);
   gbsp=GlobalBandStructDelete(gbsp);
   return gbsp;
}

/* BANDALIGN */
static GlobalBandStructPtr CC_CreatBandStruct(SeqLocPtr slp1, SeqLocPtr slp2, MashPtr msp)
{
   GlobalBandStructPtr gbsp;
   Boolean is_prot;

   is_prot = (Boolean)(msp->is_prot || msp->translate_prot);
   gbsp = CreatBandStruct(slp1, slp2, NULL,(Boolean) is_prot, (Int2) msp->band_method);
   if (( ChangeGlobalBandMatrix(gbsp, is_prot, msp->matrixname,(Int4) msp->penalty, (Int4)msp->reward)) != 0) {
      gbsp = GlobalBandStructDelete (gbsp);
      return NULL;
   }
   SetGlobaltOptions(gbsp,msp->lg1_ext,msp->rg1_ext,
		          msp->lg2_ext, msp->rg2_ext, 
   		          msp->lg1_open, msp->lg2_open, 
		          msp->rg1_open, msp->rg2_open,
		          (Int2)msp->gap_open, (Int2) msp->gap_extend);
   return gbsp;
}
/********************************************************/
static SeqAlignPtr BandAlignTwoSeqLocsFunc (SeqLocPtr slp1, SeqLocPtr slp2, MashPtr msp)
{
  GlobalBandStructPtr gbsp;
  SeqAlignPtr         newalign = NULL;
  SeqIntPtr sit;
  Boolean   is_prot;

  if (slp1==NULL || slp2==NULL)
     return NULL;
  gbsp = CC_CreatBandStruct(slp1, slp2, msp);
  if (gbsp!=NULL) 
  {
    is_prot = (Boolean)(msp->is_prot || msp->translate_prot);
    gbsp->options->low = -(Int4) SeqLocLen(slp1);
    gbsp->options->up =(Int4) SeqLocLen(slp2);
    gbsp->options->start_diag = gbsp->options->low;
    gbsp->options->width = gbsp->options->up-gbsp->options->low + 1;
/*     
    if (gbsp->options->width>2*BAND_LIMIT) {
      SetLowUpFromBlast(gbsp->options, is_prot, (Int2)0, (Int2)30,slp1, slp2);
    }
    else 
    {
       gbsp->search_type = (Uint1) G_BAND_QGAP; 
    }
*/
    newalign = GlobalBandToSeqAlign(gbsp);
    if (newalign != NULL)
    {    
       if (SeqLocStrand(slp1) == Seq_strand_minus) {
          sit=(SeqIntPtr)slp1->data.ptrvalue;
          sit->strand = Seq_strand_plus;
       }
       if (SeqLocStrand(slp2) == Seq_strand_minus) {
          sit=(SeqIntPtr)slp2->data.ptrvalue;
          sit->strand = Seq_strand_plus;
       } 
       AdjustOffSetsInSeqAlign (newalign, slp1, slp2);
    } 
    gbsp = DestructBandStruct(gbsp);
  } 
  else {
    ErrPostEx(SEV_WARNING, 0, 0, "Could not Create GlobalBandStruct");
  }
  return newalign;
}

static SeqAlignPtr BandAlignTwoSeqLocs (SeqLocPtr slp1, SeqLocPtr slp2, MashPtr msp)
{
  SeqAlignPtr         newalign = NULL;
  ValNodePtr          vnp;

  newalign = BandAlignTwoSeqLocsFunc (slp1, slp2, msp);
  if (newalign == NULL) {
     if ((msp->is_prot || msp->translate_prot)) {
        if (msp->reward > 1)
           msp->reward = 1;
        if (msp->penalty < -1)
           msp->penalty = -1;
     }
     else {
        if (msp->reward > 1)
           msp->reward = 1;
        if (msp->penalty < -1)
           msp->penalty = -1;
     }
     newalign = BandAlignTwoSeqLocsFunc (slp1, slp2, msp);
  } 
  if (newalign == NULL) {
     vnp = NULL;
     ValNodeAddPointer (&vnp, 0, slp1);
     ValNodeAddPointer (&vnp, 0, slp2);
     newalign = SeqLocToFastaSeqAlign (vnp); 
     ValNodeFree (vnp);
  }
  return newalign;
}

/********************************************************/
static Int2 number_hits (SeqAlignPtr salp)
{
  SeqAlignPtr salptmp;
  Int2 j=0;

  if (salp == NULL)
     return 0;
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next)
     j++;
  return j;
}

/*****************************************************
***
*** try_insert_seqalign
***    salplst = list of seqalign 
***    salp = tobe inserted
***    choice = 0 everything is kept
***             1 no overlaps
***             2 overlap in seq 1
***             3 overlap in seq 2
***             4 overlaps in both sequences of same strand
***             5 
***
*******************************************************/
static Boolean precede (Int4 pos1, Int4 pos2, Boolean is_plus)
{
  if (is_plus)
     return (Boolean) (pos1 < pos2);
  return (Boolean) (pos1 > pos2);
}

static SeqAlignPtr try_insert_seqalign (SeqAlignPtr salplst, SeqAlignPtr salp, Uint1 choice)
{
  SeqAlignPtr  tmp, tmp2;
  Int4         start, stop, start2, stop2,
               starti, stopi, 
               starti2, stopi2,
               startnext,
               startnext2, stopnext2;
  Boolean      goOn,
               st1, st2;

  if (salp == NULL)
     return NULL;
  if (salplst == NULL) {
     salplst = SeqAlignDup (salp);
     return salplst;
  }
  starti= SeqAlignStart (salp, 0);
  stopi = SeqAlignStop (salp, 0);
  starti2= SeqAlignStart (salp, 1);
  stopi2 = SeqAlignStop (salp, 1);
  st1 = (Boolean) !(SeqAlignStrand(salp,0) == Seq_strand_minus);
  st2 = (Boolean) !(SeqAlignStrand(salp,1) == Seq_strand_minus);

  tmp = salplst;
  start = SeqAlignStart (tmp, 0);
  stop = SeqAlignStop (tmp, 0);
  start2= SeqAlignStart (tmp, 1);
  stop2 = SeqAlignStop (tmp, 1);
  if (choice == 1)
     goOn =(Boolean) (precede(stopi, start, st1) && precede(stopi2, start2, st2));
  else if (choice == 2)
     goOn = (Boolean) (precede(stopi, start, st1) && precede(starti2, start2, st2) && precede(stopi2, stop2, st2));
  else if (choice == 3)
     goOn = (Boolean) (precede(stopi2, start2, st2) && precede(starti, start, st1) && precede(stopi, stop, st1));
  else if (choice == 4)
     goOn = (Boolean) (precede(starti, start, st1) && precede(stopi, stop, st1) && precede(starti2, start2, st2) && precede(stopi2, stop2, st2));
  else if (choice == 5)
     goOn = (Boolean) (precede(starti, start, st1));
  else
     return salplst;
  if (goOn)
  {
     tmp2 = SeqAlignDup (salp);
     tmp2->next = salplst;
     salplst = tmp2;
  }
  else {
     while (tmp != NULL)
     {
        if (choice == 1)
           goOn = (Boolean) (precede(stop, starti, st1) && precede(stop2, starti2, st2));
        else if (choice == 2)
           goOn = (Boolean) (precede(stop, starti, st1) && precede(start2, starti2, st2) && precede(stop2, stopi2, st2));
        else if (choice == 3)
           goOn = (Boolean) (precede(stop2, starti2, st2) && precede(start, starti, st1) && precede(stop, stopi, st1));
        else if (choice == 4)
           goOn = (Boolean) (precede(start, starti, st1) && precede(stop, stopi, st1) && precede(start2, starti2, st2) && precede(stop2, stopi2, st2));
        else 
           goOn=(Boolean)(precede(start, starti, st1));
        if (goOn)
        {
           if (tmp->next == NULL) {
              tmp2 = SeqAlignDup (salp);
              tmp->next = tmp2;
              break;
           }
           else {
              startnext = SeqAlignStart (tmp->next, 0);
              startnext2= SeqAlignStart (tmp->next, 1);
              stopnext2 = SeqAlignStop (tmp->next, 1);
              if (choice == 5) {
                 goOn=(Boolean)(precede(starti, startnext, st1));
              }
              else {
                 goOn=(Boolean)(precede(stopi, startnext, st1)  && precede(starti2, startnext2, st2) && precede(stopi2, stopnext2, st2)) ;
              }
              if (goOn)
              {
                 tmp2 = SeqAlignDup (salp);
                 tmp2->next = tmp->next;
                 tmp->next = tmp2;
                 break;
              }
           }
        }
        tmp = tmp->next;
        if (tmp!=NULL) {
           start = SeqAlignStart (tmp, 0);
           stop = SeqAlignStop (tmp, 0);
           start2= SeqAlignStart (tmp, 1);
           stop2 = SeqAlignStop (tmp, 1);
        }
     }
  }
  return salplst;
}

static SeqAlignPtr SortBlastHits (SeqAlignPtr salp, Int4 threshold, Uint1 choice)
{
  SeqAnnotPtr sap, sap2;
  SeqAlignPtr salpdup,
              salptmp,
              head = NULL,
              select = NULL,
              tmp = NULL, pre=NULL, preselect=NULL;
  Nlm_FloatHi bit_score;
  Nlm_FloatHi evalue;
  Int4        score, 
              number,
              totalvalmax = INT4_MAX, 
              valmax;
  Int4        gap_count = 0, 
              gap_count1= 0;
  Boolean     goOn = TRUE;
  Boolean     ok;

  if (salp == NULL || choice < 0)
     return salp;
  sap = SeqAnnotForSeqAlign (salp);
  sap2 = (SeqAnnotPtr) AsnIoMemCopy ((Pointer) sap, (AsnReadFunc) SeqAnnotAsnRead, (AsnWriteFunc) SeqAnnotAsnWrite);
  if (sap2!=NULL) {
     salpdup = sap2->data;
     sap2->data = NULL;
     SeqAnnotFree (sap2);
  }
  while (goOn) 
  { 
     valmax = 0;
     pre = select = preselect = NULL;
     for (salptmp=salpdup; salptmp!=NULL; salptmp=salptmp->next) 
     {
        GetScoreAndEvalue (salptmp, &score, &bit_score, &evalue, &number);
        if (score > valmax)
        {
           valmax = score;
           select = salptmp;
           preselect = pre;
           if (head==NULL)
              gap_count1=SeqAlignGapCount (salptmp);
        }
        else if (head==NULL && gap_count1>0 && score >= threshold) {
           gap_count=SeqAlignGapCount (salptmp);
           if (gap_count==0) {
              valmax = score;
              select = salptmp;
              preselect = pre;
              gap_count1=0;
           }
        }
        pre = salptmp;
     }
     if (valmax < threshold && head == NULL)
     {
        threshold = 20;
     }
     if (valmax >= threshold && select != NULL) 
     {
        if (head!=NULL)
           ok=(Boolean)((SeqAlignStrand(head,0) == SeqAlignStrand(select,0))
              && (SeqAlignStrand(head,1) == SeqAlignStrand(select,1)));
        else ok=TRUE;
        if (ok) {
           head = try_insert_seqalign (head, select, choice);
        }
        tmp=select;
        if (preselect==NULL) {
           salpdup = select->next;
        } else {
           preselect->next = select->next; 
        }   
        tmp->next = NULL;
        SeqAlignFree (tmp);
        totalvalmax = valmax;
     }
     else 
        goOn = FALSE;
  }
  if (salpdup!=NULL) 
     SeqAlignFree (salpdup);
  return head;
}

static SeqAlignPtr SelectBestHit (SeqAlignPtr salp, Int2 number)
{
  SeqAlignPtr salptmp,
              tmp,
              select;
  Int4        val,
              lmax=0, ltmp=0;
  Int2        j=0;

  if (salp == NULL)
     return NULL;
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next)
     j++;
  if (j==0)
     return salp;
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) {
     val = SeqAlignBestScore (salptmp);
/*
     val = SeqAlignLength (salptmp);
*/
     if (val > lmax)
        lmax = val;
  }
  if (number == 2) {
    ltmp=0;
    for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) {
       val = SeqAlignBestScore (salptmp);
/*
       val = SeqAlignLength (salptmp);
*/
       if (val < lmax && val > ltmp)
          ltmp = val;
    }
    lmax=ltmp;
  }
  j=0;
  select = NULL;
  salptmp= salp;
  while (salptmp!=NULL) {
     val = SeqAlignBestScore (salptmp);
/*
     Val = SeqAlignLength (salptmp);
*/
     if (val >= lmax) {
        if (select == NULL) {
           select = SeqAlignDup (salptmp);
           tmp = select;
        } else {
           tmp=select;
           while (tmp->next!=NULL) tmp++;
           tmp->next = SeqAlignDup (salptmp);
           tmp = tmp->next;
        }
        salptmp = salptmp->next;
        tmp->next = NULL;
        j++;
        if (j==number) break;
     }
     else {
        salptmp = salptmp->next;
/*
        tmp = salptmp;
        salptmp = salptmp->next; 
        tmp->next = NULL;
        SeqAlignFree (tmp);
*/
     } 
  } 
  select = SeqAlignExtend (select, salp);
  return select;
}

static void check_strand_inSeqalign (SeqAlignPtr salp, Uint1 strand, Int2 index)
{
  DenseSegPtr dsp;
  Int2        j; 
  Uint1Ptr    strandp;

  if (salp!=NULL) {
     dsp = (DenseSegPtr) salp->segs;
     if (dsp->strands != NULL) {
        strandp = dsp->strands;
        strandp += index;
        for (j=0; j<dsp->numseg; j++) {
           if (*strandp != strand)
              *strandp = strand;
           strandp += dsp->dim;
        }
     }   
  }
}

static SeqAlignPtr BlastTwoSeqLocs (SeqLocPtr slp1, SeqLocPtr slp2, Boolean is_prot, MashPtr msp)
{
  BLAST_OptionsBlkPtr options;
  SeqAlignPtr seqalign;
  SeqIntPtr   sit;
  Uint1       strand1, strand2;
  Boolean     delete_mash = FALSE;

  if (slp1 == NULL || slp2 == NULL)
     return NULL;
  if (msp == NULL) {
     msp = MashNew (is_prot);
     delete_mash = TRUE;
  }
  if (msp == NULL)
     return NULL;
  options = BLASTOptionNew((is_prot == FALSE) ? "blastn":"blastp",(Boolean) msp->use_gapped_blast);
  if (msp!=NULL) {
     options->gap_open= msp->gap_open;
     options->gap_extend= msp->gap_extend;
     if(is_prot) {
        if(msp->matrixname!=NULL) 
           options->matrix=StringSave(msp->matrixname);
     }
     else {
        options->reward = msp->reward;
        options->penalty= msp->penalty;
        options->wordsize = msp->wordsize;
     }
/*
     options->filter = 0;                 *//** msp->filter; **/
     options->gap_x_dropoff= msp->gap_x_dropoff;
     options->gap_x_dropoff_final= msp->gap_x_dropoff_final;
  }
  if (is_prot)
     options = BLASTOptionValidate (options, "blastp");
  else
     options = BLASTOptionValidate (options, "blastn");
  if (options == NULL)
     return NULL;
  if (!is_prot) {
     strand1=SeqLocStrand(slp1);
     strand2=SeqLocStrand(slp2);
     sit=(SeqIntPtr)slp1->data.ptrvalue;
     sit->strand = Seq_strand_both;
     sit=(SeqIntPtr)slp2->data.ptrvalue;
     sit->strand = Seq_strand_both;
  }
  seqalign = BlastTwoSequencesByLoc (slp1, slp2, NULL, options);
  if (options->wordsize==0) {
     while (seqalign==NULL && (is_prot==FALSE && options->wordsize>8)) {
        options->wordsize--;
        seqalign = BlastTwoSequencesByLoc (slp1, slp2, NULL, options);
     }
  }
  options=BLASTOptionDelete(options);
  if (delete_mash)
     MemFree (msp);
  if (!is_prot) {
     sit=(SeqIntPtr)slp1->data.ptrvalue;
     sit->strand = strand1;
     sit=(SeqIntPtr)slp2->data.ptrvalue;
     sit->strand = strand2;
  }
  return seqalign;
}


static Uint1 find_bestframe (SeqLocPtr slp, Int2 genCode)
{
  ByteStorePtr  bs;
  Int4          len;
  Int4          max = 0;
  Int4          bslen;
  Uint1         bestframe = 0, 
                frame;
  CharPtr       str;

  /* Only works for genCode == ncbieaa */
  for (frame = 1; frame <= 3; frame ++) 
  {
     bs = cds_to_pept (slp, frame, (Int2) Seq_code_ncbieaa, TRUE);
     if (bs!=NULL) {
        str = (CharPtr) BSMerge (bs, NULL);
	bslen=BSLen(bs);
        BSFree (bs);
        for (len = 0; len<bslen; len++)
           if (str[len]=='*')
              break;
        MemFree (str);
        if (len > max) {
           max = len;
           bestframe = frame;
        }
     }
  }
  return bestframe;
}

static SeqLocPtr TranslateSeqLoc (SeqLocPtr slp, Int2 genCode, Uint1 *frame)
{
  BioseqPtr     bsp = NULL;
  ByteStorePtr  bs;
  SeqLocPtr     slpp;
  Int4          length,bslen;
  
  slp = seqloc2fuzzloc (slp, TRUE, TRUE);
  *frame = find_bestframe (slp, genCode);
  bs = cds_to_pept (slp, *frame, (Int2) Seq_code_ncbieaa, TRUE);

  if(genCode != (Int2) Seq_code_ncbieaa ) {
      bslen=BSLen(bs);
      bs=BSConvertSeq(bs, (Uint1) genCode,(Uint1)Seq_code_ncbieaa,(Int4) bslen);
  }
  bsp = BioseqNew ();
  if (bsp != NULL) {
     bsp->repr = Seq_repr_raw;
     bsp->mol = Seq_mol_aa;
     bsp->seq_data_type = genCode;
     bsp->seq_data = bs;
     bsp->length = BSLen (bs);
     bsp->id = MakeNewProteinSeqId (slp, NULL);
     if (*frame==2)
        length = (Int4)((SeqLocLen(slp)-1) / (Int4)3);
     else if (*frame==3) 
        length = (Int4)((SeqLocLen(slp)-2) / (Int4)3);
     else
        length = (Int4)(SeqLocLen(slp) / (Int4)3);
     slpp = SeqLocIntNew (0, length-1, Seq_strand_plus, bsp->id);
  }
  return slpp;
}
 
/*************************************************************
*** GlobalAlignTwoSeqLocs
***   if short sequences (length < 150)
***         aligns using BandAlign only
***   else 
***         Blast the 2 seqlocs
***         if find no hit:
***            aligns using BandAlign
***         else if finds 1 hit (NOW: the longest)
***            if alignment reaches ends (3' or 5')
***               extend seqalign with the missing ends
***            else 
***               aligns unaligned seqlocs using BandAlign
***               merge seqaligns
***
*** !!!SelectBestHit : select 1 hit
***                    should select several if any
***
**************************************************************/
static SeqAlignPtr AlignExtreme5 (SeqAlignPtr salp, MashPtr msp, Int4 slpstart1, Int4 start1, Int4 slpstart2, Int4 slpstop2, Int4 start2, SeqIdPtr sip1, SeqIdPtr sip2, Uint1 strand1, Uint1 strand2)
{
  SeqAlignPtr salp2;
  SeqLocPtr   newslp1, newslp2;
  Boolean     touch_end5;
  
  if (strand2 == Seq_strand_minus) 
  {
         touch_end5 = (Boolean)((slpstart1==start1) || (slpstop2==start2));
         if (touch_end5) 
         {
            salp = SeqAlignEndExtend (salp, slpstart1, slpstop2, -1, -1, start1, start2, -1, -1, strand1, strand2);
         } 
         else 
         {
            start1 -= 1; 
            start2 += 1; 
            newslp1=SeqLocIntNew (0, start1, strand1, sip1);
            newslp2=SeqLocIntNew (start2, slpstop2, strand2, sip2);
	    msp->rg1_ext=(Int4)msp->gap_extend;
	    msp->lg1_ext=0;
	    msp->rg2_ext=(Int4)msp->gap_extend;
	    msp->lg2_ext=0;
	    msp->rg1_open=(Int4)msp->gap_open;
	    msp->lg1_open=0;
	    msp->rg2_open=(Int4)msp->gap_open;
	    msp->lg2_open=0;
            salp2 = BandAlignTwoSeqLocs (newslp1, newslp2, msp);
            salp = SeqAlignMerge (salp, salp2, FALSE);
         }
  } 
  else 
  {
         touch_end5 = (Boolean)((ABS(slpstart1-start1)<3) || (ABS(slpstart2-start2)<2));
         if (touch_end5) 
         {
            salp = SeqAlignEndExtend (salp, slpstart1, slpstart2, -1, -1, start1, start2, -1, -1, strand1, strand2);
         } 
         else 
         {
            start1 -= 1; 
            start2 -= 1; 
            newslp1=SeqLocIntNew (0, start1, strand1, sip1);
            newslp2=SeqLocIntNew (0, start2, strand2, sip2);
	    msp->rg1_ext=(Int4)msp->gap_extend;
	    msp->lg1_ext=0;
	    msp->rg2_ext=(Int4)msp->gap_extend;
	    msp->lg2_ext=0;
	    msp->rg1_open=(Int4)msp->gap_open;
	    msp->lg1_open=0;
	    msp->rg2_open=(Int4)msp->gap_open;
	    msp->lg2_open=0;
            salp2 = BandAlignTwoSeqLocs (newslp1, newslp2, msp);
            salp = SeqAlignMerge (salp, salp2, FALSE);
         }
  }
  return salp;
}

static SeqAlignPtr AlignExtreme3 (SeqAlignPtr salp, MashPtr msp, Int4 slpstop1, Int4 stop1, Int4 slpstart2, Int4 slpstop2, Int4 stop2, SeqIdPtr sip1, SeqIdPtr sip2, Uint1 strand1, Uint1 strand2)
{
  SeqAlignPtr salp2;
  SeqLocPtr   newslp1, newslp2;
  Boolean     touch_end3;
  
  if (strand2 == Seq_strand_minus) 
  {
         touch_end3 = (Boolean)((slpstop1==stop1) || (slpstart2==stop2));
         if (touch_end3) 
         {
            salp = SeqAlignEndExtend (salp, -1, -1, slpstop1, slpstart2, -1, -1, stop1, stop2, strand1, strand2);
         } 
         else 
         {
            stop1 += 1;
            stop2 -= 1;
            newslp1=SeqLocIntNew(stop1,slpstop1, strand1, sip1);
            newslp2=SeqLocIntNew(0,stop2, strand2, sip2);
	    msp->lg1_ext= (Int4)msp->gap_extend;
	    msp->rg1_ext=0;
	    msp->lg2_ext= (Int4) msp->gap_extend;
	    msp->rg2_ext=0;
	    msp->lg1_open= (Int4)msp->gap_open;
	    msp->rg1_open=0;
	    msp->lg2_open= (Int4)msp->gap_open;
	    msp->rg2_open=0;
            salp2= BandAlignTwoSeqLocs (newslp1,newslp2, msp);
            salp = SeqAlignMerge (salp, salp2, TRUE);
         } 
  } 
  else 
  {
         touch_end3 = (Boolean)((ABS(slpstop1-stop1)<3) || (ABS(slpstop2-stop2)<3));
         if (touch_end3) 
         {
            salp = SeqAlignEndExtend (salp, -1, -1, slpstop1, slpstop2, -1, -1, stop1, stop2, strand1, strand2);
         } 
         else 
         {
            stop1 += 1;
            stop2 += 1;
            newslp1=SeqLocIntNew(stop1,slpstop1, strand1, sip1);
            newslp2=SeqLocIntNew(stop2,slpstop2, strand2, sip2);
	    msp->lg1_ext=(Int4)msp->gap_extend;
	    msp->rg1_ext=0;
	    msp->lg2_ext=(Int4) msp->gap_extend;
	    msp->rg2_ext=0;
	    msp->lg1_open=(Int4)msp->gap_open;
	    msp->rg1_open=0;
	    msp->lg2_open=(Int4)msp->gap_open;
	    msp->rg2_open=0;
            salp2 = BandAlignTwoSeqLocs (newslp1, newslp2, msp);
            salp = SeqAlignMerge (salp, salp2, TRUE);
         } 
  }
  return salp;
}

static SeqAlignPtr align_extrem (SeqAlignPtr salp, SeqLocPtr slp1, SeqLocPtr slp2, MashPtr msp)
{
  SeqAlignPtr tmp;
  SeqIdPtr    sip1, sip2;
  Int4        start1, start2,
              stop1, stop2;
  Int4        slpstart1, slpstart2,
              slpstop1, slpstop2;
  Uint1       strand1, strand2;
  
  CleanStrandsSeqAlign (salp);
  sip1 = SeqLocId(slp1);
  sip2 = SeqLocId(slp2);
  strand1 = SeqAlignStrand (salp, 0);
  strand2 = SeqAlignStrand (salp, 1);
  slpstart1= SeqLocStart(slp1);
  slpstop1 = SeqLocStop(slp1);
  slpstart2= SeqLocStart(slp2);
  slpstop2 = SeqLocStop(slp2);  
  start1 = SeqAlignStart (salp, 0);
  start2 = SeqAlignStart (salp, 1);
  salp = AlignExtreme5 (salp, msp, slpstart1, start1, slpstart2, slpstop2, start2, sip1, sip2, strand1, strand2);
  check_strand_inSeqalign (salp, strand1, 0);
  check_strand_inSeqalign (salp, strand2, 1);

  if (salp!=NULL)
  {
     tmp = salp;
     while (tmp->next!=NULL)
        tmp = tmp->next;
     stop1 = SeqAlignStop (tmp, 0);
     stop2 = SeqAlignStop (tmp, 1);  
     tmp = AlignExtreme3 (tmp, msp, slpstop1,  stop1,  slpstart2, slpstop2, stop2, sip1, sip2, strand1, strand2);
     check_strand_inSeqalign (tmp, strand1, 0);
     check_strand_inSeqalign (tmp, strand2, 1);
  }
  return salp;
}

/******************************************************
static Boolean is_intron (CharPtr str, Int4 len)
{
  if (str[0]=='G' && str[1]=='T' && str[len-2]=='A' && str[len-1]=='G')
     return TRUE;
  return FALSE;
}
*********************************************************/
static FloatHi is_donor (CharPtr str, Int4 len)
{
  FloatHi one[4]={0.35, 0.35, 0.19, 0.11};
  FloatHi two[4]={0.59, 0.13, 0.14, 0.14};
  FloatHi three[4]={0.08, 0.02, 0.82, 0.08};
  FloatHi four[4]={0.01, 0.01, 1.00, 0.01};
  FloatHi five[4]={0.01, 0.01, 0.01, 1.00};
  FloatHi six[4]={0.51, 0.03, 0.43, 0.03};
  FloatHi seven[4]={0.71, 0.08, 0.12, 0.09};
  FloatHi eight[4]={0.06, 0.05, 0.84, 0.05};
  FloatHi nine[4]={0.15, 0.16, 0.17, 0.52};
  FloatHi score =1.000;
  Int4   i;
  Int4  PNTR num=NULL;

  if ((num = MemNew(len*sizeof(Int4)))==NULL) {
    /* memory error */
    return(-1);
  }

  for (i = 0; i <= len; i++){
    if (str[i]=='A')
      num[i] = 0;
    if (str[i]=='C')
      num[i] = 1;
    if (str[i]=='G')
      num[i] = 2;
    if (str[i]=='T')
      num[i] = 3;
  }
  score *= one[num[0]];
  score *= two[num[1]];
  score *= three[num[2]];
  score *= four[num[3]];
  score *= five[num[4]];
  score *= six[num[5]];
  score *= seven[num[6]];
  score *= eight[num[7]];
  score *= nine[num[8]];

  MemFree(num);
  num=NULL;

  return score;
}

static Int4 getSplicePos (CharPtr str)
{
  Int4     offset = -1;
  Int4     xcursor = 0;
  Int4     c;
  FloatHi  topscore = -FLT_MAX,
           score;
  Char     tmpstr[9];
  Int4     length;

  if (str == NULL)
    return -1;
  length = MIN(StringLen(str)/2-10, 10);
  while (xcursor <= length)
  {
      for (c = 0; c < 9; c++)
      {
        tmpstr[c] = str[xcursor+c];
      }
      if ((score=is_donor(tmpstr, 8)) > topscore)
      {
        topscore = score;
        offset = xcursor;
      }
      xcursor += 1;
  }
  if (topscore > 0.000010 && offset >=0)
    return offset+3;
  return -1;
}



static SeqAlignPtr align_btwhits (SeqAlignPtr salp, SeqIdPtr sip1, SeqIdPtr sip2, MashPtr msp)
{
  SeqAlignPtr tmp,
              newsalp,
              newsalphead, newtmp;
  SeqLocPtr   slp1, slp2;
  Int4        x1, y1, x2, y2;
  Int4        len;
  Int4        offset;
  Uint1       st1, st2;
  CharPtr     str;
  Boolean     search_intron = FALSE;

  if (salp == NULL) return NULL;
  if (msp != NULL) {
     search_intron = msp->splicing;
  }
  newsalphead = SeqAlignDup (salp);
  newtmp = newsalphead;
  tmp = salp ->next;
  while (tmp != NULL)
  {
     x1 = SeqAlignStop (newtmp, 0);
     y1 = SeqAlignStart (tmp, 0);
     x2 = SeqAlignStop (newtmp, 1);
     y2 = SeqAlignStart (tmp, 1);
     st1= SeqAlignStrand (newtmp, 0);
     st2= SeqAlignStrand (newtmp, 1);

     if (x1 + 1 == y1) {
        if (y2 == x2 + 1)
        {
           newtmp = SeqAlignMerge (newtmp, tmp, TRUE);
           tmp = tmp->next;
        }
        else if (x2 >= y2) {
           SeqAlignTrunc2 (newtmp, 0, -(abs(x2-y2+1)));
        }
        else {
           if(st1!=Seq_strand_minus) y1 -= 1; else y1+=1;
           if(st2!=Seq_strand_minus) y2 -= 1; else y2+=1;
           newtmp = SeqAlignEndExtend (newtmp, -1, -1, y1, y2, -1, -1, x1, x2, st1, st2);
           newtmp = SeqAlignMerge (newtmp, tmp, TRUE);
           tmp = tmp->next;
        }
     }
     else if (x1 >= y1) {
        SeqAlignTrunc2 (newtmp, +(abs(x1-y1+1)), 0);
     }
     else
     {
        if ((st2!=Seq_strand_minus && y2 == x2 + 1) || (st2==Seq_strand_minus && x2 == y2+1))
        {
           if (search_intron) 
           {
              str =load_seq_data(sip1, x1-5, y1+1, FALSE, &len);
              offset = getSplicePos (str);
              if (offset>= -1 && offset<len-1) 
              {
                 offset -= 6;
                 if ((offset>=-6 && offset<0) || (offset>0)) 
                 {
                    SeqAlignTrunc2 (newtmp, 0, abs(offset));
                    SeqAlignTrunc2 (tmp, abs(offset), 0);
                 }
                 x1 = SeqAlignStop (newtmp, 0);
                 y1 = SeqAlignStart (tmp, 0);
                 x2 = SeqAlignStop (newtmp, 1);
                 y2 = SeqAlignStart (tmp, 1);
                 if(st1!=Seq_strand_minus) y1 -= 1; else y1+=1;
                 if(st2!=Seq_strand_minus) y2 -= 1; else y2+=1;
                 newtmp=SeqAlignEndExtend (newtmp,-1,-1,y1, y2, -1, -1, x1, x2, st1, st2);
                 newtmp = SeqAlignMerge (newtmp, tmp, TRUE);
                 tmp = tmp->next;
              }
              MemFree (str);
           } 
           else  {
              if(st1!=Seq_strand_minus) y1 -= 1; else y1+=1; 
              if(st2!=Seq_strand_minus) y2 -= 1; else y2+=1; 
              newtmp=SeqAlignEndExtend (newtmp, -1, -1, y1, y2, -1, -1, x1, x2, st1, st2);
              newtmp = SeqAlignMerge (newtmp, tmp, TRUE);
              tmp = tmp->next;
           }
        }
        else if (st2!=Seq_strand_minus && x2 >= y2) {
            SeqAlignTrunc2 (newtmp, 0, -(abs(x2-y2+1)));
        }
        else if (st2==Seq_strand_minus && y2 >= x2) {
            SeqAlignTrunc2 (tmp, abs(y2-x2+1), 0);
        }
        else {
           slp1 = SeqLocIntNew (x1+1, y1-1, st1, sip1);
           if (st2!=Seq_strand_minus)
              slp2 = SeqLocIntNew (x2+1, y2-1, st2, sip2);
           else 
              slp2 = SeqLocIntNew (y2+1, x2-1, st2, sip2);
           newsalp = BandAlignTwoSeqLocs (slp1, slp2, msp);
           if (newsalp != NULL) 
           {
              newtmp = SeqAlignMerge (newtmp, newsalp, TRUE);
              newtmp = SeqAlignMerge (newtmp, tmp, TRUE);
           } 
           tmp = tmp->next;     
        } 
     }
  }
  return newsalphead;
}

/*************************************/
typedef struct nodehit
{
  Int2        index;
  Boolean     open;
  SeqAlignPtr salp;
  Int4        score;
  Nlm_FloatHi bit_score;
  Nlm_FloatHi evalue;
  struct nodehit PNTR     child;
  struct nodehit PNTR     next;
} NodeHit, PNTR NodeHitPtr;

static Boolean possible_child (SeqAlignPtr salp1, SeqAlignPtr salp2, Uint1 choice)
{
  Int4         start, stop, start2, stop2,
               starti, stopi, 
               starti2, stopi2;
  Boolean      goOn,
               st1, st2;

  starti= SeqAlignStart (salp1, 0);
  stopi = SeqAlignStop (salp1, 0);
  starti2= SeqAlignStart (salp1, 1);
  stopi2 = SeqAlignStop (salp1, 1);
  st1 = (Boolean) !(SeqAlignStrand(salp1,0) == Seq_strand_minus);
  st2 = (Boolean) !(SeqAlignStrand(salp1,1) == Seq_strand_minus);

  start = SeqAlignStart (salp2, 0);
  stop = SeqAlignStop (salp2, 0);
  start2= SeqAlignStart (salp2, 1);
  stop2 = SeqAlignStop (salp2, 1);
  if (choice == 1)
     goOn = (Boolean) (precede(stopi, start, st1) && precede(stopi2, start2, st2));
  else if (choice == 2)
     goOn = (Boolean) (precede(stopi, start, st1) && precede(starti2, start2, st2) && precede(stopi2, stop2, st2));
  else if (choice == 3)
     goOn = (Boolean) (precede(stopi2, start2, st2) && precede(starti, start, st1) && precede(stopi, stop, st1));
  else if (choice == 4)
     goOn = (Boolean) (precede(starti, start, st1) && precede(stopi, stop, st1) && precede(starti2, start2, st2) && precede(stopi2, stop2, st2));
  else if (choice == 5)
     goOn = (Boolean) (precede(starti, start, st1));
  else
     goOn = TRUE;
  return goOn;
}

static NodeHitPtr CreateGraph (SeqAlignPtr salp, Uint1 choice)
{
  NodeHitPtr    headnode, 
                node1, node2,
                newnode,
                child;
  SeqAlignPtr   tmp1, tmp2;
  Nlm_FloatHi   bit_score;
  Nlm_FloatHi   evalue;
  Int4          score, number;
  Int2          j;
  
  headnode = (NodeHitPtr) MemNew (sizeof(NodeHit));
  headnode->index = 1;
  headnode->salp = salp;
  headnode->child = NULL;
  headnode->next = NULL;
  node1 = headnode;
  for (tmp1=salp->next, j=2; tmp1!=NULL; tmp1=tmp1->next, j++) 
  {
     newnode = (NodeHitPtr) MemNew (sizeof(NodeHit));
     newnode->index = j;
     newnode->open = TRUE;
     newnode->salp = tmp1;
     newnode->child = NULL;
     newnode->next = NULL;
     GetScoreAndEvalue (tmp1, &score, &bit_score, &evalue, &number);
     newnode->evalue = evalue;
     newnode->score = score;
     newnode->bit_score = bit_score;
     node1->next = newnode;
     node1 = node1->next;
  }
  node1 = headnode;
  for (tmp1=salp; tmp1!=NULL; tmp1=tmp1->next) 
  {
     child = NULL;
     node2 = headnode;
     for (tmp2=salp; tmp2!=NULL; tmp2=tmp2->next) 
     {
        if (possible_child (tmp1, tmp2, choice)) {
           newnode = (NodeHitPtr) MemNew (sizeof(NodeHit));
           newnode->index = -1;
           newnode->open = TRUE;
           newnode->salp = NULL;
           newnode->child = node2;
           newnode->next = NULL;
           if (node1->child == NULL) {
              node1->child = newnode;
              child = node1->child;
           } else {
              child->next = newnode;
              child = child->next;
           }
        }
        node2 = node2->next;
     }
     node1 = node1->next;
  }
  return headnode;
}

static NodeHitPtr SimplifyGraph (NodeHitPtr headnode)
{
  NodeHitPtr node1, node2, node3, node4, node5;
  Int2       gdchild;

  for (node1 = headnode; node1!=NULL; node1 = node1->next)
  {
     for (node2=node1->child; node2!=NULL; node2 = node2->next)
     {
        node3=node2->child;
        if (node3) {
           for (node4=node3->child; node4!=NULL; node4=node4->next)
           {
              gdchild = node4->child->index;
              for (node5=node1->child; node5!=NULL; node5 = node5->next)
              {
                 if (node5->child->index == gdchild)
                 {
                    node5->open = FALSE;
                    break;
                 }
              }
           }
        }
     }
  }
  return headnode;
}


static SeqAlignPtr link_solution (ValNodePtr vnp, NodeHitPtr head, Int2 total)
{
  NodeHitPtr  nodetmp;
  SeqAlignPtr headsalp = NULL,
              tmpsalp,
              newsalp, salp;
  Int2Ptr     tab;
  Int2        j;
  Boolean     first = TRUE;

  if (vnp==NULL)
     return NULL;
  tab = (Int2Ptr) vnp->data.ptrvalue;
  if (tab==NULL)
     return NULL;
  for (j=0; j<total; j++) 
  {
     if (tab[j] > -1) 
     {
        for(nodetmp=head; nodetmp!=NULL; nodetmp=nodetmp->next) {
           if (nodetmp->index == tab[j])
              break;
        }
        if (nodetmp!=NULL && nodetmp->index == tab[j]) {
           salp=nodetmp->salp;
           if (salp!=NULL) {
              newsalp = SeqAlignDup (salp);
              if (headsalp==NULL) {
                 headsalp = newsalp;
              } else {
                 tmpsalp->next = newsalp;
              }
              tmpsalp = newsalp;
           }
        }
     }
     else break;
  }
  return headsalp;
}

static void link_solutions (ValNodePtr vnp, NodeHitPtr head, Int2 total)
{
  ValNodePtr tmp;

  for (tmp=vnp; tmp!=NULL; tmp=tmp->next)
  {
     link_solution (tmp, head, total);
  }
}

static ValNodePtr find_maxsolution (ValNodePtr vnp, NodeHitPtr head, Int2 total, float *maxscore, float *minintron, Int4 *alignlens)
{
  ValNodePtr  tmp;
  NodeHitPtr  nodetmp;
  SeqAlignPtr salp;
  Int2Ptr     tab;
  float       sum;
  Int4        start, start1, stop;
  Int4        alignlen;
  Int4        bestlen = 0;
  float       bestscore = -100.00;
  float       intronlg;
  float       bestintron = FLT_MAX;
  Int2        j;
  Boolean     first=TRUE;

  for (tmp=vnp; tmp!=NULL; tmp=tmp->next)
  {
     tab = (Int2Ptr) tmp->data.ptrvalue;
     sum = intronlg = alignlen = 0;
     for (j=0; j<total; j++) 
     {
      if (tab[j] > -1) 
      {
         for(nodetmp=head; nodetmp!=NULL; nodetmp=nodetmp->next) {
            if (nodetmp->index == tab[j])
               break;
         }
          if (nodetmp!=NULL && nodetmp->index == tab[j]) {
            salp=nodetmp->salp;
            if (salp!=NULL) {
               sum += nodetmp->score;
               start = SeqAlignStart(salp,0);
               if (first) {
                  first = FALSE;
                  start1 = start;
               }
               stop = SeqAlignStop(salp,0);
               alignlen += abs(stop - start);
/***
WriteLog ("%ld..%ld %ld %ld ", start, stop, alignlen, (long)sum);
***/
            }
         }
      }
      else {
         intronlg = (float)(abs(stop - start1) - alignlen)/(float)(j-1);
         if (sum > bestscore) {
            bestscore = sum;
         }
         if (intronlg < bestintron) {
            bestintron = intronlg;
         }
         if (alignlen > bestlen) {
            bestlen = alignlen;
         }
/***
WriteLog ("= %ld > %ld, %f > %f, %ld > %ld\n",(long)sum, (long)bestscore, intronlg, bestintron, alignlen, bestlen);
***/
         break;
      }
     }
  }
  *maxscore = bestscore;
  *minintron = bestintron;
  *alignlens = bestlen;
  return vnp;
}

static ValNodePtr get_solutions (ValNodePtr vnp, NodeHitPtr head, Int2 total, Int4 totalens)
{
  ValNodePtr  tmp, 
              bestvnp;
  NodeHitPtr  nodetmp;
  SeqAlignPtr salp;
  Int2Ptr     tab;
  float       sum;
  Int4        start, start1, stop;
  float       maxscore;
  Int4        maxalignlens, alignlens;
  Int2        index, 
              j;
  float       intronlg;
  float       minintron;
  float       bestratio=0;
  Int4        bestratio1=0;
  float       x, y, z;
  Boolean     first=TRUE;

  find_maxsolution (vnp, head, total, &maxscore, &minintron, &maxalignlens);
  bestvnp = NULL;
  index = 1;
  for (tmp=vnp; tmp!=NULL; tmp=tmp->next)
  {
    if(tmp->choice==0) 
    {
     tab = (Int2Ptr) tmp->data.ptrvalue;
     sum = intronlg = alignlens = 0;
     first = TRUE;
     for (j=0; j<total; j++) 
     {
      if (tab[j] > -1) 
      {
         for(nodetmp=head; nodetmp!=NULL; nodetmp=nodetmp->next) {
            if (nodetmp->index == tab[j])
               break;
         }
         if (nodetmp!=NULL && nodetmp->index == tab[j]) {
            salp=nodetmp->salp;
            if (salp!=NULL) {
               sum += nodetmp->score;
               start = SeqAlignStart(salp,0);
               if (first) {
                  first = FALSE;
                  start1 = start;
               }
               stop = SeqAlignStop(salp,0);
               alignlens += abs(stop - start);
            }
         }
      }
      else {
         intronlg = (float)(abs(stop - start1) - alignlens)/(float)(j-1);
         x = (float)sum / (float)maxscore;
         y = (float)intronlg / (float)minintron; 
         z = (float)alignlens / (float)maxalignlens;
         if ((Int4)(1000*x*z) > bestratio1 && (float)(x*z/y) > (float)bestratio) {
/****
WriteLog ("FFF %ld %ld =%f / %f %d = %f/  %d %d  %f    %f   %f %ld  ++ %d %d %d %d  %f  %f\n", (long)sum, (long)maxscore, x, intronlg, j-1, y, alignlens, maxalignlens, z, (float)((x*z)/y), (float)bestratio, (long)bestratio1, abs(stop - start1), alignlens, start1, stop, (float)(abs(stop - start1) - alignlens), (float)minintron );
*****/
            bestratio1 = (Int4) 1000*(x*z);
            bestratio = (float)(x*z/y);
            bestvnp = tmp;
         }
         index++;
         break;
      }
     }
    }
  }
  if (bestvnp!=NULL)
     bestvnp->choice = 5;
  return bestvnp;
}

static ValNodePtr traverseGraph (NodeHitPtr head, NodeHitPtr node, Int2 total, Int2 level, Int2Ptr tab, ValNodePtr vnp, Boolean * change)	
{
  NodeHitPtr  child,
              tmp;
  Int2Ptr     tab2;
  Int2        j;
 
  child = node->child;
  tab[level] = node->index;
  *change = TRUE;
 
  while (child!=NULL) 
  {
     tmp = child->child;
     if (child->open && tmp!=NULL) {
        vnp = traverseGraph(head, tmp, total, level+1, tab, vnp, change);
     }
     child = child->next;
  }
  if (level > 0 && *change) 
  {
     tab2 = MemNew ((size_t)((total+2)*sizeof(Int2)));
     for (j=0; j<=level; j++) {
        tab2[j] = tab[j];
     }
     for (; j<total; j++)
        tab2 [j] = -1;
     ValNodeAddPointer (&vnp, 0, (Pointer)tab2);
     tab[level] = -1;
     *change = FALSE;
  }
  return vnp;
}

static SeqAlignPtr seqalign_reverse_sorting (SeqAlignPtr salp)
{
  SeqAlignPtr tmp, tmp2,next,
              tail,
              new=NULL;

  tmp = salp;
  while (tmp->next)
  {
     tmp2=tmp; 
     next = tmp2->next;
     while (next->next) {
        tmp2=tmp2->next;
        next=next->next;
     }
     if (tmp2->next!=NULL) {
        if (new==NULL) {
           new = tmp2->next;
           tail = new;
        } else  {
           tail->next = tmp2->next;
           tail = tail->next;
        }
        tmp2->next = NULL;
     }
  }
  if (tmp)
  {
     tail->next = tmp;
  }
  return new;
}


/*************************************/
static SeqAlignPtr BlastBandAlignTwoSeqLocs (SeqLocPtr slp1, SeqLocPtr slp2, SeqAlignPtr salp, Boolean is_prot, MashPtr msp)
{
  SeqAlignPtr newsalp = NULL;
  Uint1       strand;
  Boolean     delete_mash = FALSE;

  ValNodePtr  vnp;
  NodeHitPtr  headnode, node;
  ValNodePtr  bestsol;
  Int2Ptr     tab;
  Int2        j, total;
  Boolean     change;
 
  if (salp != NULL) {
     if (msp == NULL) {
        msp = MashNew (FALSE);
        delete_mash = TRUE;
     }
     if (msp == NULL)
        return NULL;
     if (number_hits (salp) > 1) 
     {
        newsalp = SortBlastHits (salp, msp->blast_threshold, msp->choice_blastfilter);
     }
     else newsalp = salp;

     if (number_hits (newsalp) > 1) 
     {
        if (newsalp != NULL && (msp->splicing)) 
        {
           headnode = CreateGraph (newsalp, msp->choice_blastfilter);
           headnode = SimplifyGraph (headnode);
           for(total=1, node=headnode; node!=NULL; node=node->next) 
              total++; 
           tab = MemNew ((size_t)((total+2)*sizeof(Int2)));
           for(j=0; j<total; j++)
              tab[j] = -1;
           vnp = NULL;
           for (node=headnode; node!=NULL; node=node->next) {
              vnp = traverseGraph (node, node, total, 0, tab, vnp, &change);
           }
           bestsol = get_solutions (vnp, headnode, total, SeqLocLen(slp2));
           newsalp = link_solution (bestsol, headnode, total);
        }
        if (newsalp != NULL) 
        {
           if ((strand=SeqAlignStrand(newsalp, 0)) == Seq_strand_minus)
           {
              newsalp = seqalign_reverse_strand (newsalp);
              newsalp = seqalign_reverse_sorting (newsalp);
           }
           newsalp = align_btwhits (newsalp, SeqLocId(slp1), SeqLocId(slp2), msp);
        }
     }
     if (newsalp != NULL) {
        newsalp = align_extrem (newsalp, slp1, slp2, msp);
     }
     if (delete_mash)
        MemFree (msp);
  }
  return newsalp;
}

static SeqAlignPtr aaSeqAlign_to_dnaSeqAlignFunc (SeqAlignPtr PNTR salpnahead, SeqAlignPtr salpnew, ValNodePtr vnp, ValNodePtr framep)
{
  SeqAnnotPtr      sap1, sap2;
  SeqAlignPtr      salptmp;
  SeqAlignPtr      salpna;

  sap1 = SeqAnnotForSeqAlign (salpnew);
  sap2 = (SeqAnnotPtr) AsnIoMemCopy ((Pointer) sap1, (AsnReadFunc) SeqAnnotAsnRead, (AsnWriteFunc) SeqAnnotAsnWrite);
  salpna = aaSeqAlign_to_dnaSeqAlign ((SeqAlignPtr)(sap2->data), vnp, framep);
  sap1->data = NULL;
  sap2->data = NULL;
  SeqAnnotFree (sap1);
  SeqAnnotFree (sap2);
  if (salpna != NULL) {
     if (*salpnahead == NULL)
	*salpnahead = salpna;
     else {
        for (salptmp=*salpnahead; salptmp->next!=NULL; salptmp=salptmp->next)
           continue;
        salptmp->next = salpna;
     }
  }
  return *salpnahead;
}

static Boolean same_moltype (ValNodePtr list, Boolean *prot)
{
  ValNodePtr vnp;
  BioseqPtr  bsp;
  SeqLocPtr  slp;
  Boolean    is_prot;
  
  if (list==NULL)
     return FALSE;
  vnp=list;
  slp = (SeqLocPtr)vnp->data.ptrvalue;
  bsp = BioseqLockById (SeqLocId (slp));
  if (bsp!=NULL)
  {
     is_prot = ISA_aa(bsp->mol);
     BioseqUnlock (bsp);
     for (; vnp!=NULL; vnp=vnp->next)
     {
        slp = (SeqLocPtr)vnp->data.ptrvalue;
        bsp = BioseqLockById (SeqLocId (slp));
        if (bsp!=NULL) 
        {
           if (is_prot != ISA_aa(bsp->mol)) {
              ErrPostEx(SEV_ERROR, 0, 0, "One imported sequence does not have the same molecule type");
              BioseqUnlock (bsp);
              return FALSE;
           }
           BioseqUnlock (bsp);
        }
        else {
           ErrPostEx(SEV_ERROR, 0, 0, "Can not find one Bioseq");
           return FALSE;
        }
     }        
  }
  *prot = is_prot;
  return TRUE;
}

static SeqAlignPtr AlignAnyway (SeqLocPtr slp1, SeqLocPtr slp2, Boolean is_prot, MashPtr msp, Boolean message)
{
  SeqAlignPtr salp,
              salpblast;
  ValNodePtr  vnp;
  Char        st1[50], st2[50];
 
  SeqIdWrite (SeqLocId(slp1), st1, PRINTID_FASTA_SHORT, 50);
  SeqIdWrite (SeqLocId(slp2), st2, PRINTID_FASTA_SHORT, 50);
  salpblast = BlastTwoSeqLocs (slp1, slp2, is_prot, msp); 
  if (salpblast!=NULL) {
     salp = BlastBandAlignTwoSeqLocs (slp1, slp2, salpblast, is_prot, msp);
     if (salp!=NULL) {
        if (message)
           Message (MSG_OK, "%s - %s : Global alignment based on BLAST local similarity", st1, st2);
        return salp;
     }
     SeqAlignFree (salpblast);
  }
  salp = BandAlignTwoSeqLocs (slp1, slp2, msp);
  if (salp!=NULL) {
     if (message)
        Message (MSG_OK, "%s - %s : Global alignment using dynamic programing algorithm", st1, st2);
     return salp;
  }
  vnp = NULL;
  ValNodeAddPointer (&vnp, 0, slp1);
  ValNodeAddPointer (&vnp, 0, slp2);
  salp = SeqLocToFastaSeqAlign (vnp); 
  ValNodeFreeType (&vnp, TypeEmpty);
  if (salp!=NULL) {
     if (message)
        Message (MSG_OK, "Import sequence %s without alignment", st2);
     return salp;
  }
  if (message)
     ErrPostEx (SEV_WARNING, 0, 0, "No alignment");
  return NULL;
}

#ifdef SALSA_DEBUG

static void TranslatemRNA (SeqLocPtr slp1, SeqLocPtr slp2, SeqAlignPtr salp)
{
  ValNodePtr  vnp,
              tmp;
  BioseqPtr   bsp,
              pbsp;
  ByteStorePtr bs;
  SeqEntryPtr sep,
              psep,
              old;
  SeqLocPtr   slp1_2,
              slp2_2;
  SeqFeatPtr  sfp1,
              sfp2,
              sfp2_2;
  CdRegionPtr crp;
  ProtRefPtr  prp;
  Int2        code;
  Uint1       frame = 0;

  CharPtr   str;

  if (slp1==NULL || slp2==NULL)
     return;
  bsp = BioseqLockById (SeqLocId(slp2)); 
  if (bsp)
  {
     vnp = GetOrfList (bsp, 10);
     BioseqUnlock(bsp);
     if (vnp)
     {
        for (tmp=vnp; tmp!=NULL; tmp=tmp->next)
        {
           slp2_2=(SeqLocPtr)tmp->data.ptrvalue;
        }
        slp2_2=(SeqLocPtr)vnp->data.ptrvalue;
        sep = SeqMgrGetSeqEntryForData (bsp);
        if (sep) {
           code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
           if (code==0)
              code = Seq_code_ncbieaa;
           if (code) {
              crp = CreateNewCdRgn (frame, FALSE, code);
              if (crp) {
                 sfp2 = CreateNewFeature (sep, NULL, SEQFEAT_CDREGION, NULL);
                 if (sfp2) {
                    ValNodeFree (sfp2->location);
                    sfp2->location = slp2_2;
                    sfp2->data.value.ptrvalue = (Pointer)crp;
/***                    
                    bs = ProteinFromCdRegion (sfp2, TRUE); 
                    pbsp = BioseqNew ();
                    if (pbsp != NULL) 
                    {
                       pbsp->repr = Seq_repr_raw;
                       pbsp->mol = Seq_mol_aa;
                       pbsp->seq_data_type = Seq_code_ncbieaa;
                       pbsp->seq_data = bs;
                       pbsp->length = BSLen (bs);
                       bs = NULL;
                       old = SeqEntrySetScope (sep);
                       pbsp->id = MakeNewProteinSeqId (sfp2->location, NULL);
                       SeqMgrAddToBioseqIndex (pbsp);
                       SeqEntrySetScope (old);
                       psep = SeqEntryNew ();
                       if (psep != NULL) {
                         psep->choice = 1;
                         psep->data.ptrvalue = (Pointer) pbsp;
                         SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) pbsp, psep);
                       }
                       AddSeqEntryToSeqEntry (sep, psep, TRUE);
                       SetSeqFeatProduct (sfp2, pbsp);
                       prp = CreateNewProtRef ("protein", NULL, NULL, NULL);
                       if (prp) {
                          sfp2_2 = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
                          if (sfp2_2 != NULL) {
                             sfp2_2->data.value.ptrvalue = (Pointer) prp;
                          }
                       }
                    }
**/                    
/***
                    slp1_2=CopySeqLocFromSeqAlign(sfp2, SeqLocId(slp1), SeqLocId(slp2_2), salp, 0, &frame); 
***/
                 }
              }
           } 
        }
     }
  } 
}
#endif

/********************************************* 
*** SeqLocListToSeqAlign 
***    aligns the sequences from the SeqLocs list (sqloc_list) 
***    returns a SeqAlign 
***    Alignment methods:
***      1: FASTA (assume that input is FASTA aligned )
***      5: BandAlign (GlobalAlignTwoSeqLocs)
**********************************************/

extern SeqAlignPtr SeqLocListToSeqAlign (ValNodePtr sqloc_list, Int2 choice, Pointer param)
{
  ValNodePtr       vnp;
  SeqLocPtr        master_slp, 
                   slp=NULL;
  SeqLocPtr        slptmp1,
                   slptmp2;
  SeqAlignPtr      salphead= NULL,
                   salptmp = NULL,
                   salpnew = NULL,
                   salpblast = NULL,
                   salpna = NULL;
  MashPtr          msp;
  ValNodePtr       framep = NULL;
  Int2             seq_number = 0;
  Uint1            frame;
  Boolean          delete_mash = FALSE;
  Boolean          is_prot = FALSE;

  if (sqloc_list == NULL) {
      ErrPostEx(SEV_WARNING, 0, 0, "NULL SeqLocList");
     return NULL;
  }
  if (choice==2 || choice==3 || choice==4 || choice==5)
     return NULL;
  if (choice == 1) 
     return (SeqAlignPtr)SeqLocToFastaSeqAlign (sqloc_list);

  if((master_slp = (SeqLocPtr) sqloc_list->data.ptrvalue)==NULL) {
     ErrPostEx(SEV_WARNING, 0, 0, "Can not read master sequence");
     return NULL;
  }
  if (!same_moltype (sqloc_list, &is_prot)) {
     ErrPostEx(SEV_WARNING, 0, 0, "Different mol_type");
     return NULL;
  }  
  seq_number = 1;
  msp = (MashPtr)param;
  if (msp == NULL) {
     msp = MashNew (is_prot);
     delete_mash = TRUE; 
  }
  if (msp == NULL)
     return NULL;

  is_prot = (Boolean)(msp->is_prot || msp->translate_prot);
  if (msp->translate_prot) {
     slptmp1 = master_slp;
     master_slp = TranslateSeqLoc(slptmp1, (Int2) Seq_code_ncbistdaa, &frame);
     ValNodeAddInt (&framep, 1, (Int4)frame);
  }
/*   
  if(msp->is_prot && msp->matrixname==NULL) msp->matrixname="BLOSUM62"; 
*/
  for (vnp = sqloc_list->next; vnp!=NULL; vnp = vnp->next)
  {
     salpnew = NULL;
     slp = (SeqLocPtr) vnp->data.ptrvalue;
     if ( slp != NULL ) 
     {
        if (msp->translate_prot) {
           slptmp2 = slp;
           slp = TranslateSeqLoc(slptmp2, (Int2) Seq_code_ncbistdaa, &frame);
           ValNodeAddInt (&framep, 1, (Int4)frame);
        }
        if (choice == 2) {
           ObjMgrSetHold ();
/*
COMMENT    salpnew = (SeqAlignPtr) SIM2ALN (master_slp, slp, 1); */
           ObjMgrClearHold (); 
        } 
        else if (choice == 3) {
           ObjMgrSetHold ();
/*
COMMENT    salpnew=(SeqAlignPtr)SIM3ALN_choice (master_slp, slp, (Int4) 100, TRUE, TRUE); */
           ObjMgrClearHold (); 
        } 
        else if (choice == 4) {
/*
COMMENT    salpnew = (SeqAlignPtr) CSIM (master_slp, slp, 1, 0.00, NULL); */
        } 
        else if (choice == 5) {
/*
COMMENT    salpnew = SIM4ALN_choice (master_slp, slp, 1000, 8); */
        } 
        else if (choice == 6) {
           salpblast = BlastTwoSeqLocs (master_slp, slp, is_prot, msp);
           if (salpblast!=NULL) {
              salpnew = BandAlignTwoSeqLocs (master_slp, slp, msp);
           }
           SeqAlignFree (salpblast);
        } 
        else if (choice == 7) {
           salpnew = BlastTwoSeqLocs (master_slp, slp, is_prot, msp); 
           msp->choice_blastfilter = 0;
           salpnew = SortBlastHits (salpnew, msp->blast_threshold, msp->choice_blastfilter);
        } 
        else if (choice == 10) { 
           salpblast = BlastTwoSeqLocs (master_slp, slp, is_prot, msp); 
           if (salpblast!=NULL) 
           {
              salpnew = BlastBandAlignTwoSeqLocs (master_slp, slp, salpblast, is_prot, msp);
              if (salpnew==NULL)
                 SeqAlignFree (salpblast);
           }
        } 
        else if (choice == 9) {
           salpblast = BlastTwoSeqLocs (master_slp, slp, is_prot, msp); 
           if (salpblast!=NULL) {
              salpnew = AlignmRNA2genomicToSeqAlign (master_slp, slp, salpblast, msp);
              if (salpnew==NULL)
                 SeqAlignFree (salpblast);
#ifdef SALSA_DEBUG
              else
                 TranslatemRNA (master_slp, slp, salpnew);
#endif
           }
        } 
        else if (choice == 8) {
           salpblast = BlastTwoSeqLocs (master_slp, slp, is_prot, msp); 
           if (salpblast!=NULL) {
              salpnew = AlignmRNA2genomicToSeqAlign (master_slp, slp, salpblast, msp);
              if (salpnew == NULL) 
              {
                 salpnew = BlastBandAlignTwoSeqLocs (master_slp, slp, salpblast, is_prot, msp);
              }
              if (salpnew==NULL)
                 SeqAlignFree (salpblast);
           }
        }
        else if (choice == 0)
        {
           if (seq_number<3)
              salpnew = AlignAnyway (master_slp, slp, is_prot, msp, TRUE);
           else
              salpnew = AlignAnyway (master_slp, slp, is_prot, msp, FALSE);
        }
        if (salpnew != NULL) 
        {
           if (msp->translate_prot) {
              master_slp = slptmp1;
              salpna = aaSeqAlign_to_dnaSeqAlignFunc (&salpna, salpnew, sqloc_list, framep);
           }
           if (salphead == NULL) {
              salphead = salpnew;
              salptmp = salpnew;
           }
           else {
              salptmp->next = salpnew;
              salptmp = salptmp->next;
           }
           seq_number++;
        }
     } 
     else {
        ErrPostEx(SEV_WARNING, 0, 0, "Last SeqLoc was NULL");
     }
     
  }
  if (salphead != NULL) {
     if (msp->map_align)
        salphead = SeqAlignMapOnFirstSeq (salphead);
     if (seq_number > 2 && msp->multidim)
        salphead = PairSeqAlign2MultiSeqAlign (salphead);
  }
  if (salpna != NULL) {
     if (seq_number > 2 && msp->multidim)
	salpna = PairSeqAlign2MultiSeqAlign (salpna);
     if (salpna !=NULL)
        msp->transalp = salpna;
  }
  if (delete_mash)
     MemFree (msp);
  return salphead;
}

static SeqAlignPtr AlignmRNA2genomicToSeqAlign (SeqLocPtr slp1, SeqLocPtr slp2, SeqAlignPtr salpblast, MashPtr msp)
{
  SeqAlignPtr salp=NULL;
  Boolean     delete_mash = FALSE;
  
  if (salpblast != NULL) 
  {
     if (msp==NULL) {
        msp = MashNew (FALSE);
        delete_mash = TRUE;
     }
     if (msp!=NULL) {
        msp->splicing = TRUE;
        msp->choice_blastfilter = 2;
        salp = BlastBandAlignTwoSeqLocs (slp1, slp2, salpblast, FALSE, msp);
        if (delete_mash)
           MemFree (msp);
     }
  }
  return salp;
}

extern SeqLocPtr AlignmRNA2genomic (BioseqPtr bsp1, BioseqPtr bsp2)
{
  ValNodePtr  vnp = NULL;
  SeqLocPtr   slp = NULL,
              slp1 = NULL,
              slp2 = NULL;
  SeqIdPtr    sip;
  SeqAlignPtr salp,
              salpblast;
  MashPtr     msp = NULL;

  if (bsp1==NULL || bsp2==NULL)
     return NULL;
  if ((Boolean)ISA_aa(bsp1->mol) || (Boolean)ISA_aa(bsp2->mol))
     return NULL;
  sip = SeqIdFindBest (bsp1->id, 0);
  if (sip==NULL)
     return NULL;
  slp1 = SeqLocIntNew (0, bsp1->length - 1, Seq_strand_plus, sip);
  ValNodeAddPointer(&vnp, 0, (Pointer)slp1);
  sip = SeqIdFindBest (bsp2->id, 0);
  if (sip==NULL) 
     return NULL;
  slp2 = SeqLocIntNew (0, bsp2->length - 1, Seq_strand_plus, sip);
  ValNodeAddPointer(&vnp, 0, (Pointer)slp2);
  salpblast = BlastTwoSeqLocs (slp1, slp2, FALSE, msp);
  if (salpblast!=NULL)
  {
     salp = AlignmRNA2genomicToSeqAlign (slp1, slp2, salpblast, msp);
     SeqAlignFree (salpblast);
     slp = SeqLocMixFromSeqAlign (salp, NULL);
  }
  return slp;
}

