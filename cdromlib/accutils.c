/*   accutils.c
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
 * RCS $Id: accutils.c,v 6.7 1999/01/06 14:18:36 grisha Exp $
 *
 * Author:  J. Epstein
 *
 * Version Creation Date:   10/18/93
 *
 * File Description: 
 *       Utilities which make use of the Entrez "data access library"
 *
 * Modifications:  
 * --------------------------------------------------------------------------
 * $Log: accutils.c,v $
 * Revision 6.7  1999/01/06 14:18:36  grisha
 * add defines to switch ID0/ID1 usage
 *
 * Revision 6.6  1998/06/12 19:19:10  kans
 * fixed unix compiler warnings
 *
 * Revision 6.5  1998/04/28 19:29:13  shavirin
 * Fixed minor purify detected bug.
 *
 * Revision 6.4  1998/02/11 19:50:00  kans
 * FastaSeqPort takes code parameter
 *
 * Revision 6.3  1997/11/14 22:13:50  vakatov
 * [WIN32,DLL]  Added NLM_EXTERN's
 *
 * Revision 6.2  1997/09/12 15:28:26  chappey
 * Revision changes in NormalizeSeqAlign
 *
 * Revision 6.1  1997/09/04 14:14:50  chappey
 * Revision changes in NormalizeSeqAlign
 *
 * Revision 5.15  1997/08/07 16:02:40  kans
 * Revision added NormalizeSeqAlignId (Colombe)
 *
 * Revision 5.14  1997/06/26 21:55:17  vakatov
 * Revision [PC] DLL'd "ncbicdr.lib", "ncbiacc.lib", "ncbinacc.lib" and "ncbicacc.lib"
 *
 * Revision 5.13  1997/05/14 14:27:57  shavirin
 * Revision Function AccessionToFasta adopted for protein accessions
 *
 * Revision 5.11  1997/05/13 21:11:05  shavirin
 * Revision Changed function AccessionToFasta() to use PubMed accession indexes
 *
 * Revision 5.10  1997/02/04 18:58:01  epstein
 * Revision add GetFullEntrezTermList() function
 *
 * Revision 5.9  1997/01/13  15:13:34  brandon
 * changed EntrezStringToField
 *
 * Revision 5.8  1997/01/07  17:34:33  epstein
 * eliminate PMENTREZ by default
 *
 * Revision 5.7  1996/10/01  18:19:03  shavirin
 * Removed unused variable and fixed memory problems in AccessionToFasta()
 *
 * Revision 5.6  1996/09/26  17:52:45  ostell
 * made AccessionToFasta LIBCALL
 *
 * Revision 5.5  1996/09/23  21:22:09  shavirin
 * Added new function AccessionToFasta(), returning Fasta entry for
 * given Accession or GI as a string
 *
 * Revision 5.2  1996/08/14  15:15:05  brandon
 * added date parameter to tleval functions
 *
 * Revision 5.1  1996/07/01  14:06:35  epstein
 * add 'join function' EntrezCommonHierAncestor()
 *
 * Revision 4.6  1996/05/21  17:37:31  epstein
 * eliminate case-sensitivity
 *
 * Revision 4.5  1996/05/14  21:01:40  epstein
 * remove references to FLD_MLOC
 *
 * Revision 4.4  1996/04/23  19:14:46  epstein
 * memory-leak/acces cleanup, per D. Vakatov
 *
 * Revision 4.3  1995/10/11  13:40:37  epstein
 * make EntrezStringToField() and EntrezFieldToString() data-driven
 *
 * Revision 4.2  1995/10/02  02:36:00  epstein
 * add range-checking
 *
 * Revision 4.1  1995/08/24  20:44:21  epstein
 * add more stuff for genomes
 * ==========================================================================
 */

#define REVISION_STR "$Revision: 6.7 $"

#include <accutils.h>
#ifndef _CDROMLIB_
#include <cdromlib.h>
#endif

#define ERRPOST_LVL(x)               ((x) == NULL ? SEV_INFO : SEV_ERROR)

#define LEXCHAR_LPAREN                1
#define LEXCHAR_RPAREN                2
#define LEXCHAR_LBRACKET              3
#define LEXCHAR_RBRACKET              4
#define LEXCHAR_QUOTE                 5
#define LEXCHAR_AND                   6
#define LEXCHAR_OR                    7
#define LEXCHAR_NOT                   8
#define LEXCHAR_COMMA                 9
#define LEXCHAR_ATSIGN               10
#define LEXCHAR_BACKSLASH            11
#define LEXCHAR_WHITESPACE           12
#define LEXCHAR_SEMICOLON            13
#define LEXCHAR_COLON                14
#define LEXCHAR_EOL                  15
#define LEXCHAR_NULL                 16
#define LEXCHAR_OTHER                17

#define LEXSTATE_IDLE                 0
#define LEXSTATE_BACKSLASHED          1
#define LEXSTATE_INQUOTE              2
#define LEXSTATE_INQUOTE_AFTERBSLASH  3
#define LEXSTATE_INSTRING             4
#define LEXSTATE_ERROR                5

#define LEXTOK_LPAREN                 1
#define LEXTOK_RPAREN                 2
#define LEXTOK_LBRACKET               3
#define LEXTOK_RBRACKET               4
#define LEXTOK_AND                    5
#define LEXTOK_OR                     6
#define LEXTOK_NOT                    7
#define LEXTOK_COMMA                  8
#define LEXTOK_ATSIGN                 9
#define LEXTOK_STRING                10
#define LEXTOK_RANGE                 11

static Int2 lexPosition = 0;
static CharPtr lexString = NULL;
static Int2 lexState = LEXSTATE_IDLE;

static ValNodePtr nextnode = NULL;
static ValNode nextRealNode;
static Int2 lastGood = -1;
static Int2 lastBad = -1;

/****************************COLOMBE **************************/

static Int4 CCAccessionToGi (CharPtr string)
{
   CharPtr str;
   LinkSetPtr lsp;
   Int4 gi;

   str = MemNew (StringLen (string) + 10);
   sprintf (str, "\"%s\" [ACCN]", string);
   lsp = EntrezTLEvalString (str, TYP_NT, -1, NULL, NULL);
   MemFree (str);
   if (lsp == NULL) return 0;
   if (lsp->num <= 0) {
       LinkSetFree (lsp);
       return 0;
   }
   gi = lsp->uids [0];
   LinkSetFree (lsp);
   return gi;
}

/*****************************************************************************
*
*  Function:    IS_ntdb_accession
*
*  Description: Return TRUE if the input string is a validly formatted
*               nucleotide database accession number (GenBank, EMBL, DDBJ)
*
*  Arguments:   s : CharPtr; pointer to accession number string.
*                   Must be null terminated.
*
*  Author:      Mark Cavanaugh
*  Date:        7/96
*
*  WARNING:     IS_ntdb_accession() does not communicate with any central
*               resource about accession numbers. So there's no way to
*               inform it automatically about new accession number prefixes.
*
*****************************************************************************/
static Boolean IS_ntdb_accession(CharPtr s)
{
  Boolean retval = TRUE;

  Boolean first = TRUE;
  CharPtr temp;

  if (s == NULL || ! *s)
    return FALSE;

  switch (StringLen(s)) {

  case 6:                       /* Old-style 6-character accession */
    while (*s) {
      if (retval == FALSE)
        break;

      if (first) {
        if (! IS_ALPHA(*s)) {
          retval = FALSE;
          break;
        }

      switch (TO_UPPER(*s)) {
      case 'H': case 'N': case 'R': case 'T': case 'W': /* GenBank : EST */
        break;
      case 'B': case 'G': case 'I': case 'S': case 'U': /* GenBank : non-EST */
        break;
      case 'J': case 'K': case 'L': case 'M':      /* GenBank : before NCBI */
        break;
      case 'A': case 'F': case 'V': case 'X': case 'Y': case 'Z':  /* EMBL */
        break;
      case 'C': case 'D': case 'E':                /* DDBJ */
        break;
      default:
        retval = FALSE;
        break;
      }
      first = FALSE;
      }
      else {
        if (! IS_DIGIT(*s)) {
          retval = FALSE;
        }
      }
      s++;
    }  
    break;
    case 8: /* New 8-character accession, two letters + 6 digits */
 
      /* Copy the first two chars of the accession to a buffer */
 
      temp = (CharPtr) MemNew(3);

      temp[0] = *s; s++;          
      temp[1] = *s; s++;
      temp[2] = '\0';

      if ((StringICmp(temp,"AA") == 0) ||       /* NCBI EST */
          (StringICmp(temp,"AC") == 0) ||       /* NCBI HTGS */
          (StringICmp(temp,"AF") == 0) ||       /* NCBI ??? */
          (StringICmp(temp,"AD") == 0) ) {      /* NCBI accessions assigned to GSDB entries */
        /* No-op */
      }
      else if ( (StringICmp(temp,"AB") == 0) ) {        /* DDBJ */
        /* No-op */
      }
      else {
        retval = FALSE;
        break;
      }

      while (*s) {
        if (! IS_DIGIT(*s)) {
          retval = FALSE;
          break;
        }
        s++;
      }  
      break;
       
  default:   
       
    retval = FALSE;
    break;
  }                     /* Endswitch, StringLen(s) */

  return retval;                                      
}
/*****************************************************************************
*
*  Function:    IS_protdb_accession
*
*  Description: Return TRUE if the input string is a validly formatted
*               protein database accession number (SWISS-PROT)
*
*  Arguments:   s : CharPtr; pointer to accession number string.
*                   Must be null terminated.
*
*  Author:      Mark Cavanaugh
*  Date:        8/96
*
*  WARNING:     IS_protdb_accession() does not communicate with any central
*               resource about accession numbers. So there's no way to
*               inform it automatically about new accession number prefixes.
*
*****************************************************************************/
static Boolean IS_protdb_accession(CharPtr s)
{
  Boolean retval = TRUE;
  Boolean first = TRUE;

  if (s == NULL || ! *s)
    return FALSE;
 
  if (StringLen(s) != 6)
    return FALSE;

  while (*s) {
    if (retval == FALSE)
      break;

    if (first) {
      if (! IS_ALPHA(*s)) {
        retval = FALSE;
        break;
      }  

      switch (TO_UPPER(*s)) {
      case 'P': case 'Q':                 /* SWISS-PROT accessions */
        break;
      default:
        retval = FALSE;
        break;
      }  
      first = FALSE;
    }
    else {
      if (! IS_DIGIT(*s)) {
        retval = FALSE;
        break;
      }  
    }  
    s++;
  }
  return retval;
}

/****************************************************
***   Copy of Read SeqPort.bsp from SeqPort.start to stop
***   in : spp, from + to in seq coordinates
***   out: length of buffer + buffer
***        compiled with -D SAP
*****************************************************/
static Int4 AccutReadBufferFromSep (SeqPortPtr spp, CharPtr buffer, Int4 from, Int4 to, Int4 buffsegstart)
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
 

/*************************************************************
***
*** load_seq_data
***    loads bioseq sequence into a string
***    sip: SeqId of the bioseq
***    from, to: included bondaries
***    returns the length of the string (lenp)
***
***************************************************************/
static CharPtr load_seq_data (SeqIdPtr sip, Int4 from, Int4 to, Int4 *lenp)
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
     slp = SeqLocIntNew (from, to, Seq_strand_plus, sip);
     spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
     if (spp != NULL) {
        str = MemNew ((to-from+4) * sizeof(Char));
        lens = AccutReadBufferFromSep (spp, str, 0, to -from +1, 0);
        SeqPortFree (spp);
        if (lenp != NULL)
           *lenp = lens;
     }   
     SeqLocFree (slp);
  }
  return str;
}

static Int4 getlengthforid (SeqIdPtr sip)
{
  BioseqPtr        bsp;
  Int4             lens=0;
 
  if (sip==NULL)
     return 0;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     lens = bsp->length;
     BioseqUnlock (bsp);
  }
  return lens;
}

/********************************************************
***
*** NormalizeSeqAlignId
***   Checks local seqid . if a seqid string contains "acc"
***   seqid has a correspondant sequence in db.
***   This local seqid is replaced by its gi number. 
***
***   The local sequence is compared to the sequence from 
***   the db. If the local sequence is a region of the db sequence
***   the positions in the seqalign are updaded with the offset.
***
***   local copy of IS_ntdb_accession & IS_protdb_accession
***          Thanks to Mark for this useful function! 
**********************************************************/
NLM_EXTERN ValNodePtr NormalizeSeqAlignId (SeqAlignPtr salp, ValNodePtr vnp)
{
  ValNodePtr  vnptmp;
  DenseSegPtr dsp;
  SeqIdPtr    sip,
              dbsip,
              presip, 
              next,
              siptmp;
  CharPtr     TmpBuff, tmp;
  CharPtr     strlcl,
              strdb,
              strtmp;
  Char        str [52];
  Int4        gi = 0;
  Int4        offset;
  Int4        lenlcl, lendb;
  Int4        totlenlcl, totlendb;
  Int4        j, k, k1;
  Int2        index;
  Boolean     ok, goOn, found;

  Int4        step = 1000;
  Int4        window = 1020;

  EntrezInit ("Sequin", FALSE, NULL);

  if (salp!=NULL) {
     if (salp->segtype == 2) {
        dsp = (DenseSegPtr) salp->segs;
        presip = NULL;
        sip = dsp->ids;
        index = 0;
        while (sip != NULL) {
           next = sip->next;
           SeqIdWrite (sip, str, PRINTID_FASTA_LONG, 50);
           tmp = StringStr (str, "acc");
           if (tmp==NULL) {
              tmp = StringStr (str, "ACC");
           }
           if (tmp!=NULL) {
              tmp++; tmp++; tmp++;
              if (*tmp == '|')
                 tmp++;   
              TmpBuff = tmp;
              while (*tmp!='\0' && *tmp != '|' && *tmp!='\n')
                 tmp++;
              *tmp = '\0';

              ok = FALSE;
              j = StringLen (TmpBuff);
              for(k =0; k < j; k++) {
                 if(!isdigit(TmpBuff[k])) {
                    break;
                 }
              }
              if(k != j) {
                 ok=(IS_ntdb_accession(TmpBuff) || IS_protdb_accession(TmpBuff));
              }
              if (ok) {
                 gi = CCAccessionToGi(TmpBuff);
                 ok = (gi>0);
              }
              if (ok) {
                 dbsip = ValNodeNew (NULL);
                 dbsip->choice = SEQID_GI;
                 dbsip->data.intvalue = (Int4)gi;
                 sip->next = NULL;
                 found = FALSE;
                 totlenlcl = getlengthforid(sip);
                 strlcl = load_seq_data (sip, 0, MIN(totlenlcl, 20), &lenlcl); 
                 k1=0;
                 totlendb = getlengthforid(dbsip);
                 k = MIN(totlendb, k1+window); 
                 goOn = (Boolean) (strlcl!=NULL && k > 0); 
                 while (goOn) 
                 {
                    strdb = load_seq_data (dbsip, k1, k, &lendb);
                    if (strlcl != NULL && strdb != NULL) 
                    {
                     offset = 0;
                     strtmp = StringStr (strdb, strlcl);
                     if (strtmp != NULL) {
                       offset =(Int4)abs(abs((long)strdb)-abs((long)strtmp));
                       offset += k1;
                       if (offset > 0) {
                          for (j=0; j<dsp->numseg; j++) {
                             if (dsp->starts[dsp->dim*j+index] != -1)
                                dsp->starts[dsp->dim*j+index] += offset;
                          }
                       }
                       if (presip == NULL) {
                          dsp->ids = dbsip;
                       } else {
                          presip->next = dbsip;
                       }
                       dbsip->next = next;
                       vnptmp=NULL;
                       if (vnp!=NULL) {
                          for (vnptmp=vnp; vnptmp!=NULL; vnptmp=vnptmp->next) {
                             siptmp=(SeqIdPtr)vnptmp->data.ptrvalue;
                             if (SeqIdForSameBioseq(sip, siptmp))
                                break;
                          }
                       }
                       if (vnptmp==NULL)
                          ValNodeAddPointer(&vnp, 0, sip);
                       found = TRUE;
                       sip = dbsip;
                     }
                     MemFree (strdb);
                    }
                    goOn = (Boolean)(!found && k<totlendb);
                    k1=k1+step;
                    k = MIN(totlendb, k1+window);
                 }
                 if (!found)
                    sip->next = next;
                 if (strlcl != NULL)
                    MemFree (strlcl);
              } 
           }
           presip = sip;
           sip = next;
           index++;
        }
     }
  }
  EntrezFini ();
  return vnp;  
}
/****************************COLOMBE **************************/

static
Int2 LexClassifyChar(Char c)
{
    Int2 retval;

    switch(c) {
    case '(':   retval = LEXCHAR_LPAREN; break;
    case ')':   retval = LEXCHAR_RPAREN; break;
    case '[':   retval = LEXCHAR_LBRACKET; break;
    case ']':   retval = LEXCHAR_RBRACKET; break;
    case '"':   retval = LEXCHAR_QUOTE; break;
    case '&':   retval = LEXCHAR_AND; break;
    case '|':   retval = LEXCHAR_OR; break;
    case '-':   retval = LEXCHAR_NOT; break;
    case ',':   retval = LEXCHAR_COMMA; break;
    case '@':   retval = LEXCHAR_ATSIGN; break;
    case '\\':  retval = LEXCHAR_BACKSLASH; break;
    case ' ':
    case '\t':  retval = LEXCHAR_WHITESPACE; break;
    case ';':   retval = LEXCHAR_SEMICOLON; break;
    case ':':   retval = LEXCHAR_COLON; break;
    case '\0':  retval = LEXCHAR_NULL; break;
    case '\r':
    case '\n':  retval = LEXCHAR_EOL; break;
    default:    retval = LEXCHAR_OTHER; break;
    }

    return retval;
}

/* Returns -1 if no token available, else the position of the token */
static Int2
EntrezLexExpression(CharPtr str, ValNodePtr vnp)
{
    Int2 startPos;
    Int2 classChar;
    Int2 token = 0;
    Boolean done;
    Char c;
    CharPtr lexToken = NULL;
    CharPtr lexTokenStart;
    Int2 len;

    if (str == NULL && lexString == NULL)
        return -1;
    if (str != NULL)
    {
        MemFree(lexString);
        lexString = StringSave(str);
        lexPosition = 0;
        lexState = LEXSTATE_IDLE;
    }
    if (vnp == NULL)
        return -1;

    len = StringLen(lexString);
    startPos = lexPosition;

    if (lexPosition >= len)
    {
        lexState = LEXSTATE_ERROR;
        token = -1;
        lexToken = MemNew(1);
    } else {
        lexToken = MemNew(StringLen(&lexString[lexPosition]) + 1);
    }
    lexTokenStart = lexToken;

    for (done = FALSE; ! done && lexPosition <= len; lexPosition++)
    {
        c = lexString[lexPosition];
        classChar = LexClassifyChar(c);
        switch (lexState) {
        case LEXSTATE_IDLE:
            switch (classChar) {
            case LEXCHAR_LPAREN:
                token = LEXTOK_LPAREN; done = TRUE; break;
            case LEXCHAR_RPAREN:
                token = LEXTOK_RPAREN; done = TRUE; break;
            case LEXCHAR_LBRACKET:
                token = LEXTOK_LBRACKET; done = TRUE; break;
            case LEXCHAR_RBRACKET:
                token = LEXTOK_RBRACKET; done = TRUE; break;
            case LEXCHAR_AND:
                token = LEXTOK_AND; done = TRUE; break;
            case LEXCHAR_OR:
                token = LEXTOK_OR; done = TRUE; break;
            case LEXCHAR_NOT:
                token = LEXTOK_NOT; done = TRUE; break;
            case LEXCHAR_COMMA:
                token = LEXTOK_COMMA; done = TRUE; break;
            case LEXCHAR_ATSIGN:
                token = LEXTOK_ATSIGN; done = TRUE; break;
            case LEXCHAR_COLON:
                token = LEXTOK_RANGE; done = TRUE; break;
            case LEXCHAR_QUOTE:
                lexState = LEXSTATE_INQUOTE; break;
            case LEXCHAR_BACKSLASH:
                lexState = LEXSTATE_BACKSLASHED; break;
            case LEXCHAR_EOL:
            case LEXCHAR_WHITESPACE:
                startPos = lexPosition + 1; break;
            case LEXCHAR_SEMICOLON:
            case LEXCHAR_NULL:
                lexState = LEXSTATE_ERROR; done = TRUE; break;
            case LEXCHAR_OTHER:
            default:
                lexState = LEXSTATE_INSTRING; *lexToken++ = c; break;
            }
            break;
        case LEXSTATE_BACKSLASHED:
            switch (classChar) {
            case LEXCHAR_NULL:
            case LEXCHAR_EOL:
                *lexToken++ = '\0'; done = TRUE; lexState = LEXSTATE_IDLE; break;
            case LEXCHAR_LPAREN:
            case LEXCHAR_RPAREN:
            case LEXCHAR_LBRACKET:
            case LEXCHAR_RBRACKET:
            case LEXCHAR_QUOTE:
            case LEXCHAR_AND:
            case LEXCHAR_OR:
            case LEXCHAR_NOT:
            case LEXCHAR_COMMA:
            case LEXCHAR_ATSIGN:
            case LEXCHAR_BACKSLASH:
            case LEXCHAR_WHITESPACE:
            case LEXCHAR_SEMICOLON:
            case LEXCHAR_COLON:
            case LEXCHAR_OTHER:
            default:
                lexState = LEXSTATE_INSTRING; *lexToken++ = c; break;
            }
            break;
        case LEXSTATE_INQUOTE:
            switch (classChar) {
            case LEXCHAR_QUOTE:
                token = LEXTOK_STRING;
                *lexToken++ = '\0';
                done = TRUE;
                lexState = LEXSTATE_IDLE;
                break;
            case LEXCHAR_BACKSLASH:
                lexState = LEXSTATE_INQUOTE_AFTERBSLASH; break;
            case LEXCHAR_NULL:
            case LEXCHAR_EOL:
                lexState = LEXSTATE_ERROR; done = TRUE; break;
            default:
                *lexToken++ = c; break;
            }
            break;
        case LEXSTATE_INQUOTE_AFTERBSLASH:
            switch (classChar) {
            case LEXCHAR_NULL:
            case LEXCHAR_EOL:
                lexState = LEXSTATE_ERROR; done = TRUE; break;
            default:
                lexState = LEXSTATE_INQUOTE; *lexToken++ = c; break;
            }
            break;
        case LEXSTATE_INSTRING:
            switch (classChar) {
            case LEXCHAR_WHITESPACE:
            case LEXCHAR_SEMICOLON:
            case LEXCHAR_NULL:
            case LEXCHAR_EOL:
                token = LEXTOK_STRING;
                *lexToken++ = '\0';
                done = TRUE;
                lexState = LEXSTATE_IDLE;
                break;
            case LEXCHAR_BACKSLASH:
                lexState = LEXSTATE_BACKSLASHED;
                break;
            case LEXCHAR_QUOTE:
                lexState = LEXSTATE_INQUOTE;
                break;
            case LEXCHAR_OTHER:
                *lexToken++ = c; break;
            default:
                token = LEXTOK_STRING;
                *lexToken++ = '\0';
                done = TRUE;
                lexState = LEXSTATE_IDLE;
                lexPosition--; /* push back the last character */
                break;
            }
            break;
        case LEXSTATE_ERROR:
            done = TRUE;
            break;
        }
    }

    vnp->choice = (Uint1) token;
    vnp->data.ptrvalue = NULL;
    if (token == LEXTOK_STRING)
    {
        vnp->data.ptrvalue = lexTokenStart;
    } else {
        MemFree(lexTokenStart);
    }
    if (lexState == LEXSTATE_ERROR)
        return -1;
    else
        return startPos;
}

static void StrNextNode(void)
{
    nextnode = &nextRealNode;

    if (EntrezLexExpression(NULL, nextnode) < 0)
    {
        nextnode = NULL;
    }

    lastGood = lastBad;
    lastBad = lexPosition;
}

static Boolean StrExpression PROTO((ValNodePtr elst, DocType db, DocField fld));

static Boolean FindTermRequired(CharPtr term)
{
    Int2 len;

    if (term == NULL)
        return FALSE;
    len = StrLen(term);
    return ((len > 3 && term[len-1] == '.' && term[len-2] == '.' &&
        term[len-3] == '.') || (len > 1 && term[len-1] == '*') ||
        StrChr(term, '?') != NULL);
}

/* traverse all the possible fields, and take the "union" of this term */
/* over all the possibilities for the specified database               */
static Boolean
AddAllFields(ValNodePtr elst, CharPtr term, DocType db, Boolean special)
{
    EntrezInfoPtr eip;
    Boolean first = TRUE;
    Int4 specialCount;
    Int4 totalCount;
    EntrezFieldDataPtr fields;
    DocField fld;
    Boolean findTermRequired;

    if ((eip = EntrezGetInfo()) == NULL || elst == NULL)
        return FALSE;
    if (db < 0 || db >= eip->type_count)
        return FALSE;
    findTermRequired = FindTermRequired(term);
    fields = eip->types[db].fields;
    for (fld = 0; fld < eip->field_count; fld++)
    {
        if (fields == NULL || fields[fld].num_terms <= 0 || fld ==
            FLD_ORGN_HIER || (db == TYP_ML && fld == FLD_PROT))
            continue;
        if (first)
        {
            EntrezTLAddLParen(elst);
            first = FALSE;
        } else {
            EntrezTLAddOR(elst);
        }
        if (findTermRequired)
        {
            EntrezFindTerm(db, fld, term, &specialCount, &totalCount);
        }
        EntrezTLAddTerm(elst, term, db, fld, special);
    }

    if (! first)
    {
        EntrezTLAddRParen(elst);
    }

    return TRUE;
}

static Boolean
StrFactor(ValNodePtr elst, DocType db, DocField fld)
{
    if (nextnode == NULL)
    {
        ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "Null factor");
        return FALSE;
    }
    if (nextnode->choice == LEXTOK_LPAREN)
    {
        if (elst != NULL)
            EntrezTLAddLParen (elst);
        StrNextNode();
        if (! StrExpression(elst, db, fld))
        {
            return FALSE;
        }
        if (nextnode != NULL && nextnode->choice == LEXTOK_RPAREN)
        {
            if (elst != NULL)
                EntrezTLAddRParen (elst);
            StrNextNode();
        } else {
            ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "Missing right paren");
            return FALSE;
        }
    } else {
        if (nextnode->choice == LEXTOK_STRING)
        {
            CharPtr term;
            Boolean special = FALSE;
            CharPtr fldStr;
            CharPtr highRangeStr = NULL;
            Boolean allFields = FALSE;

            term = (CharPtr) (nextnode->data.ptrvalue);
            StrNextNode();
            if (nextnode != NULL && nextnode->choice == LEXTOK_RANGE)
            {
                StrNextNode();
                if (nextnode == NULL || nextnode->choice != LEXTOK_STRING)
                {
                    ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "missing second half of range");
                    MemFree(term);
                    return FALSE;
                }
                highRangeStr = (CharPtr) (nextnode->data.ptrvalue);
                StrNextNode();
            }
            if (nextnode != NULL && nextnode->choice == LEXTOK_LBRACKET)
            {
                StrNextNode();
                if (nextnode == NULL || nextnode->choice != LEXTOK_STRING)
                {
                    ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "missing field id after bracket");
                    MemFree(term);
                    MemFree(highRangeStr);
                    return FALSE;
                }
                fldStr = (CharPtr) (nextnode->data.ptrvalue);
                if (fldStr != NULL && StrCmp(fldStr, "*") == 0)
                {
                    allFields = TRUE;
                } else {
                    fld = EntrezStringToField(db, fldStr);
                }
                MemFree(nextnode->data.ptrvalue);
                if (!allFields && fld < 0)
                {
                    ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "bad field identifier");
                    MemFree(term);
                    MemFree(highRangeStr);
                    return FALSE;
                }
                StrNextNode();
                if (nextnode == NULL || (nextnode->choice != LEXTOK_COMMA &&
                    nextnode->choice != LEXTOK_RBRACKET))
                {
                    ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "missing right bracket");
                    MemFree(term);
                    MemFree(highRangeStr);
                    return FALSE;
                }
                if (nextnode->choice == LEXTOK_COMMA)
                {
                    StrNextNode();
                    if (nextnode == NULL || nextnode->choice != LEXTOK_STRING ||
                        StringCmp(nextnode->data.ptrvalue, "S") != 0)
                    {
                        ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "field qualifier error");
                        MemFree(term);
                        MemFree(highRangeStr);
                        return FALSE;
                    }
                    MemFree(nextnode->data.ptrvalue);
                    special = TRUE;
                    StrNextNode();
                    if (nextnode == NULL || nextnode->choice != LEXTOK_RBRACKET)
                    {
                        ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "missing right bracket");
                        MemFree(term);
                        MemFree(highRangeStr);
                        return FALSE;
                    }
                }
                StrNextNode();
            }

            /* if ( the default specified by caller is -1 ==> all ) then */
            if (fld < 0)
            {
                allFields = TRUE;
                if (highRangeStr != NULL)
                {
                    ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "ranges require explicit field specification");
                    MemFree(term);
                    MemFree(highRangeStr);
                    return FALSE;
                }
            }

            if (elst != NULL)
            {
                if (allFields)
                {
                   AddAllFields(elst, term, db, special);
                } else {
                    Int4 specialCount;
                    Int4 totalCount;

                    if (FindTermRequired(term))
                    {
                        EntrezFindTerm(db, fld, term, &specialCount, &totalCount);
                    }
                    EntrezTLAddTermWithRange(elst, term, db, fld, special, highRangeStr);
                }
            }
            MemFree (term);
            MemFree(highRangeStr);
        } else {
            ErrPostEx(ERRPOST_LVL(elst), ERR_CD_LEX, 0, "invalid token");
            return FALSE;
        }
    }

    return TRUE;
}


static Boolean
StrTerm(ValNodePtr elst, DocType db, DocField fld)
{
    if (! StrFactor(elst, db, fld))
        return FALSE;
    while (nextnode != NULL && nextnode->choice == LEXTOK_AND)
    {
        if (elst != NULL)
            EntrezTLAddAND(elst);
        StrNextNode();
        if (! StrFactor(elst, db, fld))
            return FALSE;
    }

    return TRUE;
}
   
static Boolean
StrDiff(ValNodePtr elst, DocType db, DocField fld)
{
    if (! StrTerm(elst, db, fld))
        return FALSE;
    while (nextnode != NULL && nextnode->choice == LEXTOK_OR)
    {
        if (elst != NULL)
            EntrezTLAddOR(elst);
        StrNextNode();
        if (! StrTerm(elst, db, fld))
            return FALSE;
    }

    return TRUE;
}
   
static Boolean
StrExpression(ValNodePtr elst, DocType db, DocField fld)
{
    if (! StrDiff(elst, db, fld))
        return FALSE;
    while (nextnode != NULL && nextnode->choice == LEXTOK_NOT)
    {
        if (elst != NULL)
            EntrezTLAddBUTNOT(elst);
        StrNextNode();
        if (! StrDiff(elst, db, fld))
            return FALSE;
    }

    return TRUE;
}

NLM_EXTERN CharPtr LIBCALL EntrezFieldToString(DocType db, DocField fld)
{
    CharPtr fldStr;
    static Char str[6];

    if (EntrezIsInited())
    {
        EntrezInfoPtr eip;

        if ((eip = EntrezGetInfo()) != NULL && fld < eip->field_count)
        {
            StrNCpy(str, eip->field_info[fld].tag, sizeof(str) - 1);
            StringUpper(str);
            return str;
        }
    }

    switch(fld)
    {
    case FLD_WORD:
        fldStr = "WORD"; break;
    case FLD_MESH:
        fldStr = "MESH"; break;
    case FLD_AUTH:
        fldStr = "AUTH"; break;
    case FLD_JOUR:
        fldStr = "JOUR"; break;
    case FLD_GENE:
        fldStr = "GENE"; break;
    case FLD_KYWD:
        fldStr = "KYWD"; break;
    case FLD_ECNO:
        fldStr = "ECNO"; break;
    case FLD_ORGN:
        fldStr = "ORGN"; break;
    case FLD_ACCN:
        fldStr = "ACCN"; break;
    case FLD_PROT:
        fldStr = "PROT"; break;
    case FLD_ORGN_HIER:
        fldStr = "HIER"; break;
    case FLD_DATE:
        fldStr = "DATE"; break;
    case FLD_FKEY:
        fldStr = "FKEY"; break;
    case FLD_PROP:
        fldStr = "PROP"; break;
    case FLD_SUBS:
        fldStr = "SUBS"; break;
    default:
        fldStr = "????";
    }

    return StringSave(fldStr);
}

NLM_EXTERN DocField LIBCALL EntrezStringToField(DocType db, CharPtr str)
{
    if (str == NULL)
        return -1;

#ifdef _PMENTREZ_
    return PMEntrezStringToField(db,str);
#else

    if (EntrezIsInited())
    {
        EntrezInfoPtr eip;
        DocField fld;

        if ((eip = EntrezGetInfo()) != NULL)
        {
            for (fld = 0; fld < eip->field_count; fld++)
            {
                if (StringICmp(str, eip->field_info[fld].tag) == 0)
                    return fld;
            }
        }


	
    }

    if (StringICmp(str, "WORD") == 0)
        return FLD_WORD;
    if (StringICmp(str, "MESH") == 0)
        return FLD_MESH;
    if (StringICmp(str, "AUTH") == 0)
        return FLD_AUTH;
    if (StringICmp(str, "JOUR") == 0)
        return FLD_JOUR;
    if (StringICmp(str, "GENE") == 0)
        return FLD_GENE;
    if (StringICmp(str, "KYWD") == 0)
        return FLD_KYWD;
    if (StringICmp(str, "ECNO") == 0)
        return FLD_ECNO;
    if (StringICmp(str, "ORGN") == 0)
        return FLD_ORGN;
    if (StringICmp(str, "ACCN") == 0)
        return FLD_ACCN;
    if (StringICmp(str, "PROT") == 0)
        return FLD_PROT;
    if (StringICmp(str, "HIER") == 0)
        return FLD_ORGN_HIER;
    if (StringICmp(str, "DATE") == 0)
        return FLD_DATE;
    if (StringICmp(str, "FKEY") == 0)
        return FLD_FKEY;
    if (StringICmp(str, "PROP") == 0)
        return FLD_PROP;
    if (StringICmp(str, "SUBS") == 0)
        return FLD_SUBS;
    return -1;
#endif
}


NLM_EXTERN LinkSetPtr LIBCALL EntrezPMTLEvalString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end, void * edc)
{
    LinkSetPtr lsp;
    ValNodePtr elst;

    if (begin != NULL)
    {
        *begin = -1;
    }
    if (end != NULL)
    {
        *end = -1;
    }

    if (str == NULL || *str == '\0')
    {
        return NULL;
    }

    if (db != TYP_ML && db != TYP_SEQ && db != TYP_AA && db != TYP_NT && db != TYP_ST && db != TYP_CH)
    {
        return NULL;
    }

    if ((elst = EntrezTLNew(db)) == NULL)
    {
        return NULL;
    }

    EntrezLexExpression(str, NULL);
    StrNextNode();

    lsp = NULL;
    if (StrExpression(elst, db, fld) && nextnode == NULL)
    {
#ifdef _PMENTREZ_
      lsp = EntrezPMTLEval(elst,edc);
#else
      lsp = EntrezTLEval(elst);
#endif /* _PMENTREZ_ */
    }
    EntrezTLFree(elst);

    if (lastGood < lastBad)
    {
        lastGood = lastBad;
    }
    if (begin != NULL)
    {
        *begin = lastGood;
    }
    if (end != NULL)
    {
        *end = lastBad;
    }

    return lsp;
}

NLM_EXTERN LinkSetPtr LIBCALL EntrezTLEvalString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end)
{
  return EntrezPMTLEvalString(str,db,fld,begin,end,NULL);
}

NLM_EXTERN ByteStorePtr LIBCALL EntrezPMTLEvalXString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end, void * edc)
{
    ByteStorePtr bsp;
    ValNodePtr elst;

    if (begin != NULL)
    {
        *begin = -1;
    }
    if (end != NULL)
    {
        *end = -1;
    }

    if (str == NULL || *str == '\0')
    {
        return NULL;
    }

    if (db != TYP_ML && db != TYP_SEQ && db != TYP_AA && db != TYP_NT && db != TYP_ST && db != TYP_CH)
    {
        return NULL;
    }

    if ((elst = EntrezTLNew(db)) == NULL)
    {
        return NULL;
    }

    EntrezLexExpression(str, NULL);
    StrNextNode();

    bsp = NULL;
    if (StrExpression(elst, db, fld) && nextnode == NULL)
    {
#ifdef _PMENTREZ_
      bsp = EntrezPMTLEvalX(elst,edc);
#else
      bsp = EntrezTLEvalX(elst);
#endif
    }

    EntrezTLFree(elst);

    if (lastGood < lastBad)
    {
        lastGood = lastBad;
    }
    if (begin != NULL)
    {
        *begin = lastGood;
    }
    if (end != NULL)
    {
        *end = lastBad;
    }

    return bsp;
}

NLM_EXTERN ByteStorePtr LIBCALL EntrezTLEvalXString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end)
{
  return EntrezPMTLEvalXString(str,db,fld,begin,end,NULL);
}

NLM_EXTERN Int4 LIBCALL EntrezPMTLEvalCountString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end, void * edc)
{
    Int4 count;
    ValNodePtr elst;

    if (begin != NULL)
    {
        *begin = -1;
    }
    if (end != NULL)
    {
        *end = -1;
    }

    if (str == NULL || *str == '\0')
    {
        return 0;
    }

    if (db != TYP_ML && db != TYP_SEQ && db != TYP_AA && db != TYP_NT && db != TYP_ST && db != TYP_CH)
    {
        return 0;
    }

    if ((elst = EntrezTLNew(db)) == NULL)
    {
        return 0;
    }

    EntrezLexExpression(str, NULL);
    StrNextNode();

    count = 0;
    if (StrExpression(elst, db, fld) && nextnode == NULL)
    {
#ifdef _PMENTREZ_
      count = EntrezPMTLEvalCount(elst,edc);
#else
      count = EntrezTLEvalCount(elst);
#endif /* _PMENTREZ_ */
    }
    EntrezTLFree(elst);

    if (lastGood < lastBad)
    {
        lastGood = lastBad;
    }
    if (begin != NULL)
    {
        *begin = lastGood;
    }
    if (end != NULL)
    {
        *end = lastBad;
    }

    return count;
}

NLM_EXTERN Int4 LIBCALL EntrezTLEvalCountString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end)
{
  return EntrezPMTLEvalCountString(str,db,fld,begin,end,NULL);
}


NLM_EXTERN Boolean LIBCALL EntrezPMTLParseString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end, void * edc)
{
    Boolean retval;

    if (begin != NULL)
    {
        *begin = -1;
    }
    if (end != NULL)
    {
        *end = -1;
    }

    if (str == NULL || *str == '\0')
    {
        return FALSE;
    }

    if (db != TYP_ML && db != TYP_SEQ && db != TYP_AA && db != TYP_NT && db != TYP_ST && db != TYP_CH)
    {
        return FALSE;
    }
    EntrezLexExpression(str, NULL);
    StrNextNode();

    retval = StrExpression(NULL, db, fld) && nextnode == NULL;

    if (lastGood < lastBad)
    {
        lastGood = lastBad;
    }
    if (begin != NULL)
    {
        *begin = lastGood;
    }
    if (end != NULL)
    {
        *end = lastBad;
    }

    return retval;
}

NLM_EXTERN Boolean LIBCALL EntrezTLParseString(CharPtr str, DocType db, DocField fld, Int2Ptr begin, Int2Ptr end)
{
  return EntrezPMTLParseString(str,db,fld,begin,end,NULL);
}


NLM_EXTERN CharPtr LIBCALL EntrezCommonHierAncestor(EntrezHierarchyPtr ehp1, EntrezHierarchyPtr ehp2, Int2Ptr distance1, Int2Ptr distance2)
{
    Int2 count;
    CharPtr retval;
 
    if (ehp1 == NULL || ehp2 == NULL || ehp1->lineage == NULL || ehp2->lineage == NULL)
        return NULL;
 
    for (count = 0; count < ehp1->numInLineage && count < ehp2->numInLineage;
         count++) {
        if (StrCmp(ehp1->lineage[count], ehp2->lineage[count]) != 0)
            break;
    }
 
    if (count <= 0)
        return NULL;
    count--;
    retval = StringSave(ehp1->lineage[count]);

    if (distance1 != NULL)
        *distance1 = ehp1->numInLineage - count;
    if (distance2 != NULL)
        *distance2 = ehp2->numInLineage - count;

    return retval;
}
 
/*****************************************************************************
*
*  Function:	AccessionToFasta
*
*  Description: Returns Fasta entry for given Accession or GI number
*
*****************************************************************************/
#define ATF_INIT_BUFF_SIZE 256

NLM_EXTERN FastaSeqPtr LIBCALL AccessionToFasta(CharPtr string) 
{
    Int4 gi=0;
    SeqEntryPtr sep;
    BioseqPtr bsp = NULL;
    SeqPortPtr spp;
    FastaSeqPtr fseq;
    Char buff[512];
    Int4 SequenceLen = 0;
    Uint1 code;
    Boolean is_na;
    
    CharPtr str;
    ByteStorePtr bstore;
    Int4 GiNum;
    
    if(string == NULL)
        return NULL;
    
    if((fseq = MemNew(sizeof(FastaSeq))) == NULL)
        return NULL;
    
    fseq->label = NULL;
    fseq->seq = NULL;
    
    if((gi = atol(string)) > 0) {
        fseq->gi = gi;
    } else {
        
        str = (CharPtr) MemNew(StringLen(string)+30);
        sprintf(str, "\"%s\"[ACCN]", string);

        if((bstore = EntrezTLEvalXString(str, TYP_NT, 
                                         -1, NULL, NULL)) == NULL ||
           (GiNum = BSLen(bstore)/sizeof(DocUid)) != 1) {
            
            /* Try protein accessions */
            
            if((bstore = EntrezTLEvalXString(str, TYP_AA, 
                                             -1, NULL, NULL)) == NULL ||
               (GiNum = BSLen(bstore)/sizeof(DocUid)) != 1) {
                
                MemFree(fseq);
                MemFree(str);
                return NULL;
            }
        }
        
        MemFree(str);
        
        BSSeek(bstore, 0L, 0);  
        
        BSRead(bstore, &gi, sizeof(Int4));
        BSFree(bstore);  
        
        fseq->gi = gi;
    }
    
    /* Now fetching sequence and defline from Entrez */
    
    if((sep = EntrezSeqEntryGet(gi, 1)) == NULL) {
        return NULL;
    }
    
    if((bsp = find_big_bioseq(sep)) == NULL) {
      SeqEntryFree(sep);
      return NULL;
    }

    StringCpy(buff, ">");
    SeqIdWrite(bsp->id, buff+1, PRINTID_FASTA_LONG, sizeof(buff));
    StringCat(buff, " ");
    CreateDefLine(NULL, bsp, buff+StringLen(buff), 
                  sizeof(buff)-StringLen(buff), 
                  0, NULL, NULL);
    
    fseq->label = StringSave(buff);  
    is_na = ISA_na (bsp->mol);
    if (is_na)
        code = Seq_code_iupacna;
    else
        code = Seq_code_ncbieaa;
    spp = FastaSeqPort(bsp, is_na, FALSE, code);
    fseq->seq = MemNew(ATF_INIT_BUFF_SIZE+1);
    
    while (FastaSeqLine(spp, fseq->seq+SequenceLen, 
                        ATF_INIT_BUFF_SIZE, TRUE)) {
        SequenceLen += ATF_INIT_BUFF_SIZE;  
        fseq->seq = Realloc(fseq->seq, SequenceLen + ATF_INIT_BUFF_SIZE + 1);
    }
    
    SeqPortFree(spp); 
    SeqEntryFree(sep);
    return fseq;
}

struct {
  CharPtr        *theMemory;
                  Int4 count;
}               state;

static          Boolean
                collectTermsProc (CharPtr term, Int4 special, Int4 total)
{
  state.theMemory[state.count++] = term;

  return TRUE;
}

NLM_EXTERN CharPtr        * LIBCALL GetFullEntrezTermList (DocType database, DocField field, Int4Ptr count)
{
  Int4            numTerms;
  Int4            numPages;
  EntrezInfoPtr   info;
  Int2 page;

  if (!EntrezIsInited ())
    return NULL;

  info = EntrezGetInfo ();
  if (info == NULL || info->type_count <= database)
    return NULL;

  numTerms = info->types[database].fields[field].num_terms;
  numPages = info->types[database].fields[field].num_bucket;
  state.theMemory = (CharPtr *) MemNew (sizeof (CharPtr) * numTerms);
  state.count = 0;

  for (page = 0; page < numPages; page += INT2_MAX / 2) {
    EntrezTermListByPage (database, field, page, MIN ((numPages - page), INT2_MAX / 2), collectTermsProc);
  }

  *count = state.count;
  return state.theMemory;
}
