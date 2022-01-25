/* ===========================================================================
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
*  any work or product based on this material.
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
* File Name:  spideymain.c
*
* Author:  Sarah Wheelan
*
* Version Creation Date:   5/01
*
* $Revision: 6.1 $
*
* File Description: main functions for running Spidey as a standalone 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: spideymain.c,v $
* Revision 6.1  2001/05/24 16:27:41  wheelan
* initial checkin
*
*
* ==========================================================================
*/


#include <ncbi.h>
#include <spidey.h>
#include <accid1.h>
#include <lsqfetch.h>

#define MYARGGENFILE   0
#define MYARGMRNAFILE  1
#define MYARGOUTFILE   2
#define MYARGPRALIGN   3
#define MYARGALNFILE   4
#define MYARGGILIST    5
#define MYARGNUMMOD    6
#define MYARGORG       7
#define MYARG1STEVAL   8
#define MYARG2NDEVAL   9
#define MYARG3RDEVAL   10
#define MYARGIDCUT     11
#define MYARGLENCUT    12
#define MYARGSPEC      13
#define MYARGASN       14
#define MYARGASNFILE   15
/*#define MYARGDRAFTFILE 16*/
#define MYARGMASKED    16
#define MYARGGETCDS    17
#define MYARGTABLE     18
/*#define MYARGACEDB     19*/
#define MYARGFROM      19 
#define MYARGTO        20

#define NUMARGS        21

Args myargs[NUMARGS] = {
   {"Input file -- genomic sequence(s)", NULL, NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
   {"Input file -- mRNA sequence(s)", NULL, NULL, NULL, FALSE, 'm', ARG_FILE_IN, 0.0, 0, NULL},
   {"Output file", "stdout", NULL, NULL, TRUE, 'o', ARG_STRING, 0.0, 0, NULL},
   {"Print alignment?", "F", NULL, NULL, TRUE, 'p', ARG_BOOLEAN, 0.0, 0, NULL},
   {"Alignment file", "spidey.aln", NULL, NULL, TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
   {"Input file is a GI list", "F", NULL, NULL, TRUE, 'G', ARG_BOOLEAN, 0.0, 0, NULL},
   {"Number of gene models", "1", NULL, NULL, TRUE, 'n', ARG_INT, 0.0, 0, NULL},
   {"Organism (genomic sequence) v=vertebrate,\nd = drosophila, p = plant, c = C. elegans", "v", NULL, NULL, TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
   {"First-pass e-value", "0.0000001", NULL, NULL, TRUE, 'e', ARG_FLOAT, 0.0, 0, NULL},
   {"Second-pass e-value", "0.001", NULL, NULL, TRUE, 'f', ARG_FLOAT, 0.0, 0, NULL},
   {"Third-pass e-value", "10", NULL, NULL, TRUE, 'g', ARG_FLOAT, 0.0, 0, NULL},
   {"% identity cutoff", "0", NULL, NULL, TRUE, 'c', ARG_INT, 0.0, 0, NULL},
   {"% length coverage cutoff", "0", NULL, NULL, TRUE, 'l', ARG_INT, 0.0, 0, NULL},
   {"interspecies alignment", "F", NULL, NULL, TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
   {"Print ASN.1 alignment?", "F", NULL, NULL, TRUE, 'j', ARG_BOOLEAN, 0.0, 0, NULL},
   {"File for asn.1", "spidey.asn", NULL, NULL, TRUE, 'k', ARG_STRING, 0.0, 0, NULL},
/*   {"File of draft boundary information", NULL, NULL, NULL, TRUE, 'b', ARG_FILE_IN, 0.0, 0, NULL},*/
   {"Is the input mRNA masked (lowercase)?", "F", NULL, NULL, TRUE, 'w', ARG_BOOLEAN, 0.0, 0, NULL},
   {"Fetch the CDS and compute its results also?", "F", NULL, NULL, TRUE, 'd', ARG_BOOLEAN, 0.0, 0, NULL},
   {"File with feature table", NULL, NULL, NULL, TRUE, 't', ARG_FILE_IN, 0.0, 0, NULL},
   /*{"ACEDB format", "F", NULL, NULL, TRUE, 'A', ARG_BOOLEAN, 0.0, 0, NULL},*/
   {"Start of genomic interval desired (from)", "0", NULL, NULL, TRUE, 'F', ARG_INT, 0.0, 0, NULL},
   {"Stop of genomic interval desired (to)", "0", NULL, NULL, TRUE, 'T', ARG_INT, 0.0, 0, NULL},
};

static void SPI_FindAllNuc(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);
static CharPtr ReadALine (CharPtr str, size_t size, FILE *fp);
static BioseqPtr SPI_GetBspFromGIOrAcc(CharPtr str);
static void SPI_ReadFeatureTable(FILE *ifp, SPI_bsinfoPtr spim_head);

Int2 Main()
{
   AsnIoPtr           aip;
   BioseqPtr          bsp;
   Pointer            dataptr;
   Uint2              datatype;
   Boolean            found;
   SPI_mRNAToHerdPtr  h_head;
   SPI_mRNAToHerdPtr  h_prev;
   SPI_mRNAToHerdPtr  hptr;
   FILE               *ifp;
   Boolean            isGIlist;
   Char               line[60];
   Boolean            lowercase;
   SeqLocPtr          lcaseloc;
   FILE               *ofp;
   FILE               *ofp2;
   SeqAlignPtr        sap;
   SeqAnnotPtr        sanp;
   SeqEntryPtr        sep;
   SeqIdPtr           sip;
   SeqLocPtr          slp;
   SPI_bsinfoPtr      spig;
   SPI_bsinfoPtr      spig_head;
   SPI_bsinfoPtr      spig_prev;
   SPI_bsinfoPtr      spim;
   SPI_bsinfoPtr      spim_head;
   SPI_bsinfoPtr      spim_prev;
   SPI_OptionsPtr     spot;
   CharPtr            str;
   CharPtr            txt;

   ID1BioseqFetchEnable("spidey", FALSE);
   LocalSeqFetchInit(FALSE);
   /* standard setup */
   ErrSetFatalLevel (SEV_MAX);
   ErrClearOptFlags (EO_SHOW_USERSTR);
   UseLocalAsnloadDataAndErrMsg ();
   ErrPathReset ();
   if (! AllObjLoad ())
   {
      Message (MSG_FATAL, "AllObjLoad failed");
      return 1;
   }
   if (! SubmitAsnLoad ())
   {
      Message (MSG_FATAL, "SubmitAsnLoad failed");
      return 1;
   }
   if (! FeatDefSetLoad ())
   {
      Message (MSG_FATAL, "FeatDefSetLoad failed");
      return 1;
   }
   if (! SeqCodeSetLoad ())
   {
      Message (MSG_FATAL, "SeqCodeSetLoad failed");
      return 1;
   }
   if (! GeneticCodeTableLoad ())
   {
      Message (MSG_FATAL, "GeneticCodeTableLoad failed");
      return 1;
   }
   if (!GetArgs("SPIDEY", NUMARGS, myargs))
      return 0;
   /* set the error message level high to suppress warnings from BLAST */
   isGIlist = (Boolean)myargs[MYARGGILIST].intvalue;
   txt = myargs[MYARGGENFILE].strvalue;
   ifp = FileOpen(txt, "r");
   spig_head = NULL;
   if (ifp == NULL)
   {
      bsp = SPI_GetBspFromGIOrAcc(txt);
      if (bsp == NULL)
      {
         ErrPostEx(SEV_ERROR, 0, 0, "Can't open genomic input file\n");
         return -1;
      } else
      {
         spig_head = (SPI_bsinfoPtr)MemNew(sizeof(SPI_bsinfo));
         spig_head->bsp = bsp;
      }
   }
   if (spig_head == NULL)
   {
      spig_prev = NULL;
      /* read in the genomic sequence(s) first and put them into bsinfo structures */
      while ((dataptr = ReadAsnFastaOrFlatFile (ifp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL)
      {
         if (datatype == OBJ_BIOSEQ)
         {
            spig = (SPI_bsinfoPtr)MemNew(sizeof(SPI_bsinfo));
            spig->bsp = (BioseqPtr)dataptr;
            if (spig_head == NULL)
               spig_head = spig_prev = spig;
            else
            {
               spig_prev->next = spig;
               spig_prev = spig;
            }
         } else if (datatype == OBJ_SEQENTRY)
         {
            sep = (SeqEntryPtr)dataptr;
            SeqEntryExplore(sep, &spig_head, SPI_FindAllNuc);
         }
      }
      FileClose(ifp);
   }
   if (spig_head == NULL)
   {
      ErrPostEx(SEV_ERROR, 0, 0, "No valid bioseqs in genomic file\n");
      return -1;
   } else if (ISA_aa(spig_head->bsp->mol))
   {
      ErrPostEx(SEV_ERROR, 0, 0, "At least one of the genomic sequences appears to be a protein.\n");
      return -1;
   }
   if (spig_head->next != NULL)
   {
      ErrPostEx(SEV_ERROR, 0, 0, "This version can only process one genomic sequence at a time.  Only the first sequence in this file will be used.\n");
      spig_head->next = NULL;
   }
   spim_head = spim_prev = NULL;
   txt = myargs[MYARGMRNAFILE].strvalue;
   ifp = FileOpen(txt, "r");
   if (ifp == NULL)
   {
      bsp = SPI_GetBspFromGIOrAcc(txt);
      if (bsp == NULL)
      {
         ErrPostEx(SEV_ERROR, 0, 0, "Can't open mRNA input file\n");
         return -1;
      } else
      {
         spim_head = (SPI_bsinfoPtr)MemNew(sizeof(SPI_bsinfo));
         spim_head->bsp = bsp;
      }
   }
   if (spim_head == NULL)
   {
      lowercase = (Boolean)myargs[MYARGMASKED].intvalue;
      lcaseloc = NULL;
      /* if the mRNA has lowercase masking, read it in carefully to record the masking */
      if (lowercase == TRUE)
      {
         while ((sep = FastaToSeqEntryForDb(ifp, TRUE, NULL, TRUE, NULL, NULL, &lcaseloc)) != NULL)
         {
            SeqEntryExplore(sep, &spim_head, SPI_FindAllNuc);
            if (lcaseloc != NULL)  /* put masking info into the bsinfo structure */
            {
               spim = spim_head;
               sip = SeqLocId(lcaseloc);
               found = FALSE;
               while (spim != NULL && !found)
               {
                  if (SeqIdComp(sip, spim->bsp->id) == SIC_YES)
                  {
                     found = TRUE;
                     spim->lcaseloc = lcaseloc;
                  }
                  spim = spim->next;
               }
               lcaseloc = NULL;
            }
         }
      } else if (isGIlist) /* mRNA file is a list of GIs, must fetch the bioseqs */
      {
         str = ReadALine(line, sizeof(line), ifp);
         while (str != NULL)
         {
            bsp = SPI_GetBspFromGIOrAcc(str);
            if (bsp != NULL)
            {
               spim = (SPI_bsinfoPtr)MemNew(sizeof(SPI_bsinfo));
               spim->bsp = bsp;
               if (spim_head == NULL)
                  spim_head = spim_prev = spim;
               else
               {
                  spim_prev->next = spim;
                  spim_prev = spim;
               }
            }
            str = ReadALine(line, sizeof(line), ifp);
         }
      } else /* mRNAs are FASTA or ASN.1, read them all in */
      {
         while ((dataptr = ReadAsnFastaOrFlatFile (ifp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL)
         {
            if (datatype == OBJ_BIOSEQ)
            {
               spim = (SPI_bsinfoPtr)MemNew(sizeof(SPI_bsinfo));
               spim->bsp = (BioseqPtr)dataptr;
               if (spim_head == NULL)
                  spim_head = spim_prev = spim;
               else
               {
                  spim_prev->next = spim;
                  spim_prev = spim;
               }
            } else if (datatype == OBJ_SEQENTRY)
            {
               sep = (SeqEntryPtr)dataptr;
               SeqEntryExplore(sep, &spim_head, SPI_FindAllNuc);
            }
         }
      }
      FileClose(ifp);
   }
   if (spim_head == NULL)
   {
      ErrPostEx(SEV_ERROR, 0, 0, "No valid bioseqs in mRNA file\n");
      return -1;
   } else if (ISA_aa(spim_head->bsp->mol))
   {
      ErrPostEx(SEV_ERROR, 0, 0, "At least one of the mRNA sequences appears to be a protein\n");
      return -1;
   }
   txt = myargs[MYARGTABLE].strvalue;
   if (txt != NULL)
   {
      ifp = FileOpen(txt, "r");
      if (ifp == NULL)
      {
         ErrPostEx(SEV_ERROR, 0, 0, "Unable to open table file\n");
         return -1;
      }
      SPI_ReadFeatureTable(ifp, spim_head);
      spim = spim_head;
      while (spim != NULL)
      {
         if (spim->lcaseloc != NULL)
         {
            slp = (SeqLocPtr)ValNodeNew(NULL);
            slp->choice = SEQLOC_MIX;
            slp->data.ptrvalue = (Pointer)spim->lcaseloc;
            spim->lcaseloc = slp;
         }
         spim = spim->next;
      }
   }
   spim = spim_head;
   txt = myargs[MYARGOUTFILE].strvalue;
   ofp = FileOpen(txt, "w");
   if (ofp == NULL)
   {
      ErrPostEx(SEV_ERROR, 0, 0, "Unable to open output file\n");
      return -1;
   }
   ErrSetMessageLevel(SEV_MAX);
   spot = (SPI_OptionsPtr)MemNew(sizeof(SPI_Options));
   spot->firstpasseval = myargs[MYARG1STEVAL].floatvalue;
   spot->secpasseval = myargs[MYARG2NDEVAL].floatvalue;
   spot->thirdpasseval = myargs[MYARG3RDEVAL].floatvalue;
   spot->numreturns = myargs[MYARGNUMMOD].intvalue;
   spot->idcutoff = myargs[MYARGIDCUT].intvalue;
   spot->lencutoff = myargs[MYARGLENCUT].intvalue;
   spot->printaln = (Boolean)myargs[MYARGPRALIGN].intvalue;
   spot->interspecies = (Boolean)myargs[MYARGSPEC].intvalue;
   spot->printasn = (Boolean)myargs[MYARGASN].intvalue;
   spot->fetchcds = (Boolean)myargs[MYARGGETCDS].intvalue;
   /*spot->ace = (Boolean)myargs[MYARGACEDB].intvalue;*/
   spot->from = myargs[MYARGFROM].intvalue;
   spot->to = myargs[MYARGTO].intvalue;
   txt = myargs[MYARGORG].strvalue;
   if (!StringICmp(txt, "d") || !StringICmp(txt, "D"))
      spot->organism = SPI_FLY;
   else if (!StringICmp(txt, "p") || !StringICmp(txt, "P"))
      spot->organism = SPI_PLANT;
   else if (!StringICmp(txt, "c") || !StringICmp(txt, "C"))
      spot->organism = SPI_CELEGANS;
   else
      spot->organism = SPI_VERTEBRATE;
   if (spot->printaln)
   {
      txt = myargs[MYARGALNFILE].strvalue;
      ofp2 = FileOpen(txt, "a");
   } else
      ofp2 = NULL;
   sap = NULL;
   if (spot->printasn)
      spot->sap_head = &sap;
   /*txt = myargs[MYARGDRAFTFILE].strvalue;
   if (txt != NULL)
      spot->draftfile = StringSave(txt);*/
   h_head = h_prev = NULL;
   while (spim != NULL)
   {
      spot->lcaseloc = spim->lcaseloc;
      if (spot->draftfile == NULL)
         SPI_AlnSinglemRNAToGen(spig_head, spim, ofp, ofp2, spot);
      else
      {
         hptr = SPI_AlnSinglemRNAToPieces(spig_head, spim, ofp, ofp2, spot);
         if (h_head != NULL)
         {
            h_prev->next = hptr;
            h_prev = hptr;
         } else
            h_head = h_prev = hptr;
      }
      spim = spim->next;
   }
   /* create the ASN.1 output, if requested; need to use the continuous alignment */
   /* that was generated */
   if (spot->printasn && *(spot->sap_head) != NULL && spot->draftfile == NULL)
   {
      sanp = SeqAnnotForSeqAlign(*(spot->sap_head));
      txt = myargs[MYARGASNFILE].strvalue;
      aip = AsnIoOpen(txt, "w");
      SeqAnnotAsnWrite(sanp, aip, NULL);
      AsnIoClose(aip);
      SeqAlignSetFree(*(spot->sap_head));
   }
   FileClose(ofp);
   FileClose(ofp2);
   SPI_OptionsFree(spot);
   SPI_bsinfoFreeList(spim_head);
   SPI_bsinfoFreeList(spig_head);
   LocalSeqFetchDisable();
   ID1BioseqFetchDisable();
   return 0;
}

/***************************************************************************
*
*  SPI_FindAllNuc is a SeqEntryExplore callback function that puts all
*  the nucleotide bioseqs in a seqentry into a linked list of SPI_bsinfoPtrs.
*
***************************************************************************/
static void SPI_FindAllNuc(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
   BioseqPtr      bsp;
   SPI_bsinfoPtr  spi;
   SPI_bsinfoPtr  PNTR spim_head;

   spim_head = (SPI_bsinfoPtr PNTR)data;
   spi = *spim_head;
   if (IS_Bioseq(sep))
   {
      bsp = (BioseqPtr)sep->data.ptrvalue;
      if (ISA_na(bsp->mol))
      {
         if (spi != NULL)
         {
            while (spi->next != NULL)
            {
               spi = spi->next;
            }
            spi->next = (SPI_bsinfoPtr)MemNew(sizeof(SPI_bsinfo));
            spi = spi->next;
         } else
            spi = *spim_head = (SPI_bsinfoPtr)MemNew(sizeof(SPI_bsinfo));
         spi->bsp = bsp;
      }
   }
}

/***************************************************************************
*
*  ReadALine is copied from Jonathan's code -- it simply reads 
*  consecutive lines from a file.
*
***************************************************************************/
static CharPtr ReadALine (CharPtr str, size_t size, FILE *fp)
{
   Char     ch;
   CharPtr  ptr;
   CharPtr  rsult;

   if (str == NULL || size < 1 || fp == NULL) return NULL;
   *str = '\0';
   rsult = fgets (str, size, fp);
   if (rsult != NULL)
   {
      ptr = str;
      ch = *ptr;
      while (ch != '\0' && ch != '\n' && ch != '\r')
      {
         ptr++;
         ch = *ptr;
      }
      *ptr = '\0';
   }
   return rsult;
}

/***************************************************************************
*
*  SPI_GetBspFromGIOrAcc takes a string, decides whether it's probably
*  a GI or an accession, converts it to an appropriate seqid, and then
*  attempts to fetch and return the bioseq referred to by that id.
*
***************************************************************************/
static BioseqPtr SPI_GetBspFromGIOrAcc(CharPtr str)
{
   BioseqPtr  bsp;
   Int4       gi;
   Char       ptr;
   SeqIdPtr   sip;
   ValNode    vn;

   str = TrimSpacesAroundString(str);
   ptr = *str;
   if (IS_ALPHA(ptr))  /* accession */
   {
      sip = SeqIdFromAccessionDotVersion(str);
      bsp = BioseqLockById(sip);
   } else  /* it's a GI */
   {
      gi = atoi(str);
      vn.choice = SEQID_GI;
      vn.data.intvalue = gi;
      vn.next = NULL;
      bsp = BioseqLockById(&vn);
   }
   return bsp;
}

/***************************************************************************
*
*  SPI_ReadFeatureTable reads in a tab-delimited file and converts it
*  to feature information:
*
*  sequence id	name of feature		start	stop
*  NM_004377.1	repetitive_region	12	40
*
*  Masking and other feature information for all mRNAs should be in the
*  same file; ids not in the mRNA list or unknown feature names will be
*  ignored.
*
***************************************************************************/
static void SPI_ReadFeatureTable(FILE *ifp, SPI_bsinfoPtr spim_head)
{
   Char           c;
   CharPtr        fields[5];
   Boolean        found;
   Char           line[60];
   Int4           i;
   CharPtr        ptr;
   SeqIdPtr       sip;
   SeqLocPtr      slp;
   SPI_bsinfoPtr  spim;
   CharPtr        str;
   ValNode        vn;

   str = ReadALine(line, sizeof(line), ifp);
   while (str != NULL)
   {
      ptr = strtok(str, " \t");
      for (i=0; i<4; i++)
      {
         fields[i] = StringSave(ptr);
         ptr = strtok(NULL, " \t");
      }
      if (fields[1] != NULL && !StrCmp(fields[1], "repetitive_region"))
      {
         c = *fields[0];
         if (IS_ALPHA(c))
            sip = SeqIdFromAccessionDotVersion(fields[0]);
         else
         {
            vn.choice = SEQID_GI;
            vn.data.intvalue = atoi(fields[0]);
            vn.next = NULL;
            sip = &vn;
         }
         spim = spim_head;
         found = FALSE;
         while (!found && spim != NULL)
         {
            if (SeqIdIn(sip, spim->bsp->id) == TRUE)
               found = TRUE;
            else
               spim = spim->next;
         }
         if (found)
         {
            slp = SeqLocIntNew(atoi(fields[2]), atoi(fields[3]), Seq_strand_both, spim->bsp->id);
            slp->next = spim->lcaseloc;
            spim->lcaseloc = slp;
         }
      }
      str = ReadALine(line, sizeof(line), ifp);
   }
}
