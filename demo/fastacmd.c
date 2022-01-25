/* $Id: fastacmd.c,v 6.7 1998/02/11 18:06:43 madden Exp $
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
* File Name:  $RCSfile: fastacmd.c,v $
*
* Author:  Sergei Shavirin
*
* Initial Version Creation Date: 05/20/1997
*
* $Revision: 6.7 $
*
* File Description:
*        FASTA retrievel system using ISAM indexes
*
* $Log: fastacmd.c,v $
* Revision 6.7  1998/02/11 18:06:43  madden
* Fix for reading IDs in from a file
*
* Revision 6.6  1998/02/06 18:26:35  madden
* Removed stripping of white spaces
*
* Revision 6.5  1998/01/29 19:47:02  madden
* Changed second call from BioseqRawToFasta to BioseqRawToFastaExtra
*
* Revision 6.4  1998/01/27 20:27:03  madden
* Added option to specify sequence line length
*
* Revision 6.3  1998/01/23 21:56:16  madden
* Error messages sent to stderr
*
* Revision 6.2  1998/01/16 22:04:20  madden
* Call to readdb_new_ex, fixed FUM
*
* Revision 6.1  1997/11/07 16:19:15  shavirin
* Added possibility to retrieve redundant accessions
*
* Revision 6.0  1997/08/25 18:19:56  madden
* Revision changed to 6.0
*
* Revision 5.2  1997/05/20 21:00:45  shavirin
* Remove spurious error message in accession/locus retrieval
*
* Revision 5.1  1997/05/20 15:47:30  shavirin
* Initial revision
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <readdb.h>


#define BUFFER_LENGTH 128

typedef struct FCMDAccList {
    CharPtr acc;
    Int4 gi;
    struct FCMDAccList *next;
} FCMDAccList, PNTR FCMDAccListPtr;

#define NUMARG 5

static Args myargs [NUMARG] = {
    { "Database", 
      "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
    { "Search string: GIs, accessions and locuses may be used delimited\n"
      "      by comma or space)",
      NULL, NULL, NULL, TRUE, 's', ARG_STRING, 0.0, 0, NULL},
    { "Input file wilth GIs/accessions/locuses for batch retrieval",
      NULL, NULL, NULL, TRUE, 'i', ARG_STRING, 0.0, 0, NULL},
    { "Retrieve duplicated accessions",
      "F", NULL, NULL, TRUE, 'a', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Line length for sequence", 
      "80", NULL, NULL, TRUE, 'l', ARG_INT, 0.0, 0, NULL}
};

static FCMDAccListPtr GetAccList(CharPtr file, Int4Ptr TotalItems);
static Boolean IS_protdb_accession(CharPtr s);
static Boolean IS_ntdb_accession(CharPtr s);
static void FCMDAccListFree(FCMDAccListPtr falp);

Int2 Main (void)
 
{
  BioseqPtr bsp;
  ReadDBFILEPtr rdfp;
  Int4 i, fid, TotalItems=0, count;
  FCMDAccListPtr falp, falp_tmp;
  CharPtr buffer = NULL;
  FILE *fd;
  Int4Ptr ids = NULL;

  if (! GetArgs ("fastacmd", NUMARG, myargs)) {
      return (1);
  }

  if( !ErrSetLogfile ("stderr", ELOG_APPEND) ) {
      exit(1);
  }

  if((rdfp = readdb_new_ex(myargs [0].strvalue, READDB_DB_UNKNOWN, FALSE)) == NULL) {
      fprintf(stderr, "ERROR: Cannot initialize readdb engine. Exiting...\n");
      return(1);
  }

  if(myargs[1].strvalue != NULL) {
      if((falp =  GetAccList(myargs[1].strvalue, &TotalItems)) == NULL) {
          fprintf(stderr, "ERROR: No valid Gis/Accessions found. Exiting...\n");
          return(1);
      }
  } else if(myargs[2].strvalue != NULL){
      if((fd = FileOpen(myargs[2].strvalue, "r")) == NULL)
          return (1);

      buffer = WWWReadFileInMemory(fd, 0, TRUE);

      if((falp =  GetAccList(buffer, &TotalItems)) == NULL) {
          fprintf(stderr, "ERROR: No valid Gis/Accessions found. Exiting...\n");
          return(1);
      }
  }

  /*  printf("**** %d Gis/Accessions present in "
      "Fastacmd request\n\n", TotalItems); */
  
  for (falp_tmp = falp; falp_tmp != NULL; falp_tmp = falp_tmp->next) {  
      
      if(falp_tmp->gi != 0) {
          fid = readdb_gi2seq(rdfp, falp_tmp->gi);
      } else {
          if(myargs[3].intvalue == 0) {
              fid = readdb_acc2fasta(rdfp, falp_tmp->acc);
          } else {
              count = 0;
              fid = readdb_acc2fastaEx(rdfp, falp_tmp->acc, &ids, &count);
          }
      }
      
      if(fid < 0 && fid != -1) { 
          fprintf(stderr, "Accesion search failed for \"%s\" "
                 "with error code %d\n", falp_tmp->acc, fid);
          return 1;
      } else if (fid == -1) {
          fprintf(stderr, "Entry \"%s\" not found\n", falp_tmp->acc);
      } else if (ids == NULL) { /* gi or SeqId */
          bsp = readdb_get_bioseq(rdfp, fid);   
          BioseqRawToFastaExtra(bsp, stdout, myargs[4].intvalue);
          BioseqFree(bsp);  
      } else {
          for(i = 0; i < count; i++) {
              bsp = readdb_get_bioseq(rdfp, ids[i]);   
              BioseqRawToFastaExtra(bsp, stdout, myargs[4].intvalue);
              BioseqFree(bsp);
          }
      }   
  }
  readdb_destruct(rdfp);
  MemFree(buffer);
  FCMDAccListFree(falp);
  return 0;
} 

static FCMDAccListPtr GetAccList(CharPtr file, 
                                 Int4Ptr TotalItems)
{
    Char TmpBuff[128];
    Int4 i, j, k;
    Int4 FileLen = 0;
    FCMDAccListPtr AccList = NULL;
    FCMDAccListPtr AccListTmp, AccListLast;
    Int4 NumNotValid = 0;
    Int4 gi = 0;
  
  if(file == NULL || file[0] == NULLB) {
    *TotalItems = 0;
    return NULL;
  }
  
  FileLen = StringLen(file);
  
  for(i = 0; i < FileLen; i++) {
      
      if(isspace(file[i]) || file[i] == ',') /* Rolling spaces */
          continue;
      
      /* This is defence from badly formatted requests */
      
      if(NumNotValid > 10) {
          fprintf(stderr, "**** ERROR: Too many invalid Gis/Accessions, "
                 "parsing aborted\n");
          *TotalItems = 0;
          return NULL;
      }
      
      /* Rolling spaces */
      
      j= 0;
      while (j < 128  && i < FileLen) { 
          TmpBuff[j] = file[i];
          j++; i++;
          if(file[i] == ',')  /* Comma is valid delimiter */
              break;
      }
      TmpBuff[j] = NULLB;
    
      /* Is gi/accession too long ??? */

      if(j == 128) {
          fprintf(stderr, "**** WARNING: Gi/Accession \"%s\" is too long\r\n", 
                 TmpBuff);
          NumNotValid++;
          
          while(!isspace(file[i]) || 
                file[i] == ',' || 
                file[i] == NULLB) /* Rolling until spaces */
              i++;
          continue;  /* Next may be valid ... who knows...?? */   
      }
      
      /* Now validating accession/gi */
      
      for(k =0; k < j; k++) {
          if(!IS_DIGIT(TmpBuff[k])) {
              break;
          }
      }

      if(k == j)
          gi = atol(TmpBuff);
      
      /* If this is valid Accession check and tranfer it to gi */
      
      /* It we come here - we got valid text ID */
      
      if(AccList == NULL) { /* first element */
          AccList = (FCMDAccListPtr) MemNew(sizeof(FCMDAccList));
          AccListTmp = AccList;
          AccListTmp->acc = StringSave(TmpBuff);
          AccListTmp->gi = gi;
          AccListTmp->next = NULL;
          AccListLast=AccListTmp;
          *TotalItems = *TotalItems +1; 
      } else {
          AccListTmp = (FCMDAccListPtr) MemNew(sizeof(FCMDAccList));
          AccListLast->next = AccListTmp;
          AccListTmp->acc = StringSave(TmpBuff);
          AccListTmp->gi = gi;
          AccListTmp->next = NULL;
          AccListLast = AccListTmp;
          *TotalItems = *TotalItems +1;
      }
  }
  if(NumNotValid) {
      fprintf(stderr, "**** %d invalid Gi%s/Accession%s present in fastacmd "
             "request\r\n", 
             NumNotValid,
             NumNotValid == 1 ? "" : "s",
             NumNotValid == 1 ? "" : "s"
             );
  }
  return AccList;
}

static void FCMDAccListFree(FCMDAccListPtr falp)
{
    FCMDAccListPtr falp_tmp, falp_next;

    if(falp == NULL)
        return;

    for(falp_tmp = falp; falp_tmp != NULL; falp_tmp=falp_next) {
        falp_next = falp_tmp->next;
        MemFree(falp_tmp->acc);
        MemFree(falp_tmp);
    }
}

/*****************************************************************************
*
*  Function:	IS_ntdb_accession
*
*  Description:	Return TRUE if the input string is a validly formatted
*		nucleotide database accession number (GenBank, EMBL, DDBJ)
*
*  Arguments:	s : CharPtr; pointer to accession number string.
*		    Must be null terminated.
*
*  Author:	Mark Cavanaugh
*  Date:	7/96
*
*  WARNING:	IS_ntdb_accession() does not communicate with any central
*		resource about accession numbers. So there's no way to
*		inform it automatically about new accession number prefixes.
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
    
  case 6:			/* Old-style 6-character accession */
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
      case 'J': case 'K': case 'L': case 'M':	   /* GenBank : before NCBI */
        break;
      case 'A': case 'F': case 'V': case 'X': case 'Y': case 'Z':  /* EMBL */
        break;
      case 'C': case 'D': case 'E':		   /* DDBJ */
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
      
      if ((StringICmp(temp,"AA") == 0) ||	/* NCBI EST */
	  (StringICmp(temp,"AC") == 0) ||	/* NCBI HTGS */
	  (StringICmp(temp,"AD") == 0) ) {	/* NCBI accessions assigned to GSDB entries */
        /* No-op */
      }
      else if ( (StringICmp(temp,"AB") == 0) ) {	/* DDBJ */
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
  }			/* Endswitch, StringLen(s) */
  
  return retval;
}
/*****************************************************************************
*
*  Function:	IS_protdb_accession
*
*  Description:	Return TRUE if the input string is a validly formatted
*		protein database accession number (SWISS-PROT)
*
*  Arguments:	s : CharPtr; pointer to accession number string.
*		    Must be null terminated.
*
*  Author:	Mark Cavanaugh
*  Date:	8/96
*
*  WARNING:	IS_protdb_accession() does not communicate with any central
*		resource about accession numbers. So there's no way to
*		inform it automatically about new accession number prefixes.
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
      case 'P': case 'Q':		  /* SWISS-PROT accessions */
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

