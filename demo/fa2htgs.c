/*****************************************************************************
*
*   fa2htgs.c
*      This is a framework for reading an ASN.1 Seq-entry or Bioseq-set,
*       doing some processing on it, and outputting it again. seqport.h
*       is included since it covers most necessary utilities that you might
*       need. You may need to add others for specialized reports and so on.
*
*       The check for aipout == NULL is to show how to change the code if
*         no output is desired. Change the default in myargs from "stdout" to
*         NULL to make output command line optional.
*
*       This program can be used "as is" to convert between binary and text
*        ASN.1 through the object loaders.
*
*   -- original version from Jim Ostell
*   -- version 1.3, add "-m" and "-u" parameters   Hsiu-Chuan  2-27-97
*   -- version 1.4, add "-x" parameter             Hsiu-Chuan  4-10-97 
*   -- test ci to CVS in ncbi demo directory       Hsiu-Chuan  5-1-97
*   -- version 1.5, add "-C" for clone-lib
*                       "-M" for map
*                       "-O" for comment from a file
*                       "-T" for phrap file          Hsiu-Chuan  1-30-98
*                       "-P" for contig names
*                       "-A" for accession list file  Kans  4-8-98
*                       "-X" for coordinates on resulting master or individual accessions
*
*****************************************************************************/
#include <subutil.h>
#include <tofasta.h>
#include <sqnutils.h>

#define NUMARG 26
Args myargs[NUMARG] = {
   {"Filename for fasta input","stdin",NULL,NULL,TRUE,'i',ARG_FILE_IN,0.0,0,NULL},
   {"Filename for Seq-submit template","template.sub",NULL,NULL,FALSE,'t',ARG_FILE_IN,0.0,0,NULL},
   {"Filename for asn.1 output","stdout", NULL,NULL,TRUE,'o',ARG_FILE_OUT,0.0,0,NULL},
   {"Log errors to file named:",NULL,NULL,NULL,TRUE,'e',ARG_FILE_OUT,0.0,0,NULL} ,
   {"Organism name?","Homo sapiens", NULL ,NULL ,TRUE,'n',ARG_STRING,0.0,0,NULL},
   {"Sequence name?",NULL, NULL ,NULL ,FALSE,'s',ARG_STRING,0.0,0,NULL},
   {"length of sequence in bp?","0", NULL ,NULL ,FALSE,'l',ARG_INT,0.0,0,NULL},
   {"Genome Center tag?",NULL, NULL ,NULL ,FALSE,'g',ARG_STRING,0.0,0,NULL},
   {"HTGS phase?","1", "0" ,"3" ,FALSE,'p',ARG_INT,0.0,0,NULL},
   {"GenBank accession (if an update)",NULL, NULL ,NULL ,TRUE,'a',ARG_STRING,0.0,0,NULL},
   {"Remark for update?",NULL, NULL ,NULL ,TRUE,'r',ARG_STRING,0.0,0,NULL},
   {"Clone name?",NULL, NULL ,NULL ,TRUE,'c',ARG_STRING,0.0,0,NULL},
   {"Chromosome?",NULL, NULL ,NULL ,TRUE,'h',ARG_STRING,0.0,0,NULL},
   {"Title for sequence?",NULL, NULL ,NULL ,TRUE,'d',ARG_STRING,0.0,0,NULL},
   {"Take comment from template ?","F", NULL ,NULL ,TRUE,'m',ARG_BOOLEAN,0.0,0,NULL},
   {"Take biosource from template ?","F", NULL ,NULL ,TRUE,'u',ARG_BOOLEAN,0.0,0,NULL},
   {"Secondary accession number, separate by comas if multiple, s.t. U10000,L11000", NULL, NULL ,NULL ,TRUE,'x',ARG_STRING,0.0,0,NULL},
   {"Clone library name?",NULL, NULL ,NULL ,TRUE,'C',ARG_STRING,0.0,0,NULL},
   {"Map?",NULL, NULL ,NULL ,TRUE,'M',ARG_STRING,0.0,0,NULL},
   {"Filename for the comment:",NULL,NULL,NULL,TRUE,'O',ARG_FILE_IN,0.0,0,NULL} ,
   {"Filename for phrap input",NULL,NULL,NULL,TRUE,'T',ARG_FILE_IN,0.0,0,NULL} ,
   {"Contigs to use, separate by comas if multiple", NULL, NULL ,NULL ,TRUE,'P',ARG_STRING,0.0,0,NULL},
   {"Filename for accession list input",NULL,NULL,NULL,TRUE,'A',ARG_FILE_IN,0.0,0,NULL} ,
   {"Coordinates are on the resulting sequence ?","F", NULL ,NULL ,TRUE,'X',ARG_BOOLEAN,0.0,0,NULL},
   {"HTGS_DRAFT sequence?","F", NULL ,NULL ,TRUE,'D',ARG_BOOLEAN,0.0,0,NULL},
   {"Strain name?",NULL, NULL ,NULL ,TRUE,'S',ARG_STRING,0.0,0,NULL},
};

/*------------- MakeAc2GBSeqId() -----------------------*/
/***************************************************************
*   MakeAc2GBSeqId:
*   -- return NULL if acnum == null
*                                             Hsiu-Chuan 4-18-97
****************************************************************/
static SeqIdPtr  MakeAc2GBSeqId(CharPtr accession)
{
   TextSeqIdPtr tsip;
   SeqIdPtr sip;

   if (accession == NULL || *accession == '\0')
      return NULL;

   sip = ValNodeNew(NULL);
   sip->choice = SEQID_GENBANK;
   tsip = TextSeqIdNew();
   sip->data.ptrvalue = tsip;
   tsip->accession = StringSave(accession);

   return sip;

} /* MakeAc2GBSeqId */

/*----------- AddExtraAc2Entry() ----------------------------*/
/***************************************************************
*   AddExtraAc2Entry:
*                                             Hsiu-Chuan 4-11-97
****************************************************************/
static Boolean AddExtraAc2Entry (SeqEntryPtr entry , CharPtr extra_ac )
{
   BioseqPtr  bsp;
   ValNodePtr vnp;
   GBBlockPtr gbp;
   Char       acnum[17];
   CharPtr    p;
   Int4       i, j;
   SeqHistPtr shp;
   SeqIdPtr   sip;

   if ((entry == NULL) || (extra_ac == NULL))
      return FALSE;

   bsp = (BioseqPtr)(entry->data.ptrvalue);

   for (gbp= NULL, vnp = bsp->descr; vnp != NULL; vnp = vnp->next)
   {
       if (vnp->choice == Seq_descr_genbank)
       {
          gbp = vnp->data.ptrvalue;
          break;
       }
   }

   shp = bsp->hist; 

   if (gbp == NULL)
   {
      vnp = (ValNodePtr) NewDescrOnSeqEntry (entry, Seq_descr_genbank);
      gbp = GBBlockNew();
      vnp->data.ptrvalue = (Pointer)gbp;
   }
   
   for (p = extra_ac; *p != '\0';)
   {
       for (i = 0; isalnum(*p) && *p != '\0'; ++p, ++i)
           acnum[i] = *p;
       acnum[i] = '\0'; 
               /* check one_letter+5digits or two_letter+6digits */
       if (i == 6 || i == 8)
       {
          if (!isalpha(acnum[0]) || (!(isdigit(acnum[1]) && i == 6) &&
              !(isalpha(acnum[1]) && i == 8)))
          {
             ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
             return FALSE;
          }

          for (j = 2; j < i; ++j)
          {
              if (!(isdigit(acnum[j])))
              {
                 ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
                 return FALSE;
              }
          }

          ValNodeCopyStr(&gbp->extra_accessions, 0, acnum);
          sip = MakeAc2GBSeqId (acnum);
          if (shp == NULL)
          {
             shp = SeqHistNew();
             bsp->hist = shp;
          }
          ValNodeLink(&shp->replace_ids, sip);
       }
       else
       {
          ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
          return FALSE;
       }

       while (!isalnum(*p) && *p != '\0')
           ++p;
   }

   return TRUE;

} /* AddExtraAc2Entry */

static void RescueSeqGraphs (BioseqPtr bsp, Int2 index, ValNodePtr PNTR vnpp)

{
  SeqAnnotPtr   nextsap;
  SeqGraphPtr   nextsgp;
  Pointer PNTR  prevsap;
  Pointer PNTR  prevsgp;
  SeqAnnotPtr   sap;
  SeqGraphPtr   sgp;

  if (bsp == NULL || vnpp == NULL) return;
  sap = bsp->annot;
  prevsap = (Pointer PNTR) &(bsp->annot);
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      prevsgp = (Pointer PNTR) &(sap->data);
      while (sgp != NULL) {
        nextsgp = sgp->next;
        *(prevsgp) = sgp->next;
        sgp->next = NULL;
        ValNodeAddPointer (vnpp, index, (Pointer) sgp);
        sgp = nextsgp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static Boolean HasNoText (CharPtr str)

{
  Char  ch;

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static SeqAnnotPtr NewSeqAnnotType3 (CharPtr name, SeqGraphPtr sgp)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (! HasNoText (name)) {
    ValNodeAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}

static void OffsetAndLinkSeqGraph (BioseqPtr bsp, SeqGraphPtr sgp, Int2 index)

{
  DeltaSeqPtr  dsp;
  SeqGraphPtr  lastsgp;
  Int4         len;
  SeqLitPtr    litp;
  SeqAnnotPtr  sap;
  SeqIntPtr    sintp;
  SeqLocPtr    slp;

  if (bsp == NULL || sgp == NULL || index < 1) return;
  len = 0;
  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext);
         dsp != NULL && index > 1; dsp = dsp->next, index--) {
      if (dsp->choice == 1) {
        len += SeqLocLen ((SeqLocPtr) dsp->data.ptrvalue);
      } else if (dsp->choice == 2) {
        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          len += litp->length;
        }
      }
    }
  }
  slp = sgp->loc;
  if (slp != NULL && slp->choice == SEQLOC_INT) {
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp != NULL) {
      sintp->from += len;
      sintp->to += len;
      sintp->id = SeqIdFree (sintp->id);
      sintp->id = SeqIdDup (bsp->id);
    }
  }
  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 3) {
      for (lastsgp = sap->data; lastsgp->next != NULL; lastsgp = lastsgp->next) {
        continue;
      }
      lastsgp->next = sgp;
      break;
    }
  }
  if (sap == NULL) {
    if (bsp->annot != NULL) {
      for (sap = bsp->annot; sap->next != NULL; sap = sap->next) {
        continue;
      }
      sap->next = NewSeqAnnotType3 ("Graphs", sgp);
    } else {
      bsp->annot = NewSeqAnnotType3 ("Graphs", sgp);
    }
  }
}

/*****************************************************************************
*
*   Main program loop to read, process, write SeqEntrys
*
*****************************************************************************/
Int2 Main(void)
{
   AsnIoPtr aip;
   FILE *fp, *cfp;
   SeqSubmitPtr ssp;
   NCBISubPtr nsp;
   SeqEntryPtr sep, oldsep, the_entry, sep_list, nextsep;
   BioseqPtr bsp, the_bsp;
   Uint1 htgs_phase; /* a value from 0-3 */
   Uint1 MI_htgs_phase;  /* mapping of htgs_phase to MI_TECH_htgs_? */
   CharPtr  newstr, accession, remark, center, organism, clone, seqbuf,
      seqname, strain, chromosome, title, extra_ac, clone_lib, map,
      comment_fname, comment_fstr, phrap_fname, fasta_fname, contigs, accn_fname;
   Char  instr[120];
   Int4   totalen, filelen, len, length = 0, cumlength = 0;
   SeqLitPtr slp;
   Int2 errs;
   BioseqSetPtr bssp;
   ValNodePtr vnp, PNTR prevpnt, next;
   Boolean   temp_org, temp_comment, lastwasraw, coordsOnMaster, htgsDraft;
   Int2 index = 0;
   ValNodePtr rescuedsgps = NULL;

   CharPtr tool_ver = "fa2htgs 1.7";

               /* check command line arguments */

   if ( ! GetArgs(tool_ver, NUMARG, myargs))
      return 1;

   fasta_fname = myargs[0].strvalue;
   organism = myargs[4].strvalue;
   seqname = myargs[5].strvalue;
   length = myargs[6].intvalue;
   center = myargs[7].strvalue;
   htgs_phase = (Uint1)(myargs[8].intvalue);
   if (htgs_phase == 0)
     MI_htgs_phase = (Uint1)MI_TECH_htgs_0;
   else
     MI_htgs_phase = (Uint1)(MI_TECH_htgs_1 + htgs_phase - 1);

   accession = myargs[9].strvalue;
   remark = myargs[10].strvalue;
   clone = myargs[11].strvalue;
   chromosome = myargs[12].strvalue;
   title = myargs[13].strvalue;
   temp_comment = (Boolean) myargs[14].intvalue;
   temp_org = (Boolean) myargs[15].intvalue;
   extra_ac = myargs[16].strvalue;
   clone_lib = myargs[17].strvalue;
   map = myargs[18].strvalue;
   comment_fname = myargs[19].strvalue;
   comment_fstr = NULL;
   phrap_fname = myargs[20].strvalue;
   contigs = myargs[21].strvalue;
   accn_fname = myargs[22].strvalue;
   coordsOnMaster = (Boolean) myargs[23].intvalue;
   htgsDraft = (Boolean) myargs[24].intvalue;
   strain = myargs[25].strvalue;

   UseLocalAsnloadDataAndErrMsg (); /* finds data directory without a .ncbirc file */

               /* load the sequence alphabets  */
               /* (and sequence parse trees)   */
   if (! SeqEntryLoad())
   {
      ErrShow();
      return 1;
   }
   if (! SubmitAsnLoad())
   {
      ErrShow();
      return 1;
   }
                                /* log errors instead of die */
   if (myargs[3].strvalue != NULL)
   {
      if (! ErrSetLog (myargs[3].strvalue))
      {
         ErrShow();
         return 1;
      }
      else
         ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
   }
                         /* open phrap contig file */
   if (phrap_fname != NULL)
   {
      if ((fp = FileOpen (phrap_fname, "r")) == NULL)
      {
         ErrPostEx(SEV_ERROR,0,0, "Can't open %s", phrap_fname);
         ErrShow();
         return 1;
      }
   }
 
   else if (accn_fname != NULL) {
      if ((fp = FileOpen (accn_fname, "r")) == NULL)
      {
         ErrPostEx(SEV_ERROR,0,0, "Can't open %s", accn_fname);
         ErrShow();
         return 1;
      }
   }
                         /* open input fasta file */
   else if ((fp = FileOpen (fasta_fname, "r")) == NULL)
   {
      ErrPostEx(SEV_ERROR,0,0, "Can't open %s", fasta_fname);
      ErrShow();
      return 1;
   }
                        /* open comment file */
   if (comment_fname != NULL) {

      if ((cfp = FileOpen (comment_fname, "r")) == NULL)
      {
         ErrPostEx(SEV_ERROR,0,0, "Can't open %s", comment_fname);
         ErrShow();
         return 1;
      }
      /* rules for building the comment string a file:
         -- maximum 100 characters per line
         -- insert a " " to concatnate lines
         -- insert a "~" to concatnate lines if
            it is a blank line or there are leading
            spaces in the beginning of the line
                                 Hsiu-Chuan Chen 1-30-98
      */ 

      filelen = FileLength (comment_fname);
      filelen = filelen + 1000;
      comment_fstr = MemNew (filelen);
      totalen = 0;
      while (fgets (instr, 110, cfp) != NULL)
      {
          len = StringLen (instr);
          while (len > 0 && instr[len-1] == '\n')
          {
             instr[len-1] = '\0';
             len = StringLen (instr);
          }

          totalen = totalen + len + 2;
          if (totalen > filelen)
          {
             filelen = filelen + 1000;
             newstr = MemNew (filelen);
             StringCpy (newstr, comment_fstr);
             MemFree (comment_fstr);

             comment_fstr = newstr;
          }

          if (comment_fstr != NULL)
          {
             if (instr[0] == '\0' || instr[0] == ' ')
                StringCat (comment_fstr, "~");
             else
                StringCat (comment_fstr, " ");

             StringCat (comment_fstr, instr);
          }
          else
             StringCpy (comment_fstr, instr);

      } /* while */
   }
                        /* open template file */
   if ((aip = AsnIoOpen (myargs[1].strvalue, "r")) == NULL)
   {
      ErrPostEx(SEV_ERROR,0,0, "Can't open %s", myargs[1].strvalue);
      ErrShow();
      FileClose(fp);
      return 1;
   }
   ssp = SeqSubmitAsnRead(aip, NULL);
   AsnIoClose(aip);
   if (ssp == NULL)
   {
      ErrPostEx(SEV_ERROR,0,0, "Can't read %s", myargs[1].strvalue);
      ErrShow();
      FileClose(fp);
      return 1;
   }

   oldsep = (SeqEntryPtr)(ssp->data);  /* clear out template */
   ssp->data = NULL;
   MemFree(ssp->sub->tool);
   ssp->sub->tool = StringSave(tool_ver);
   nsp = MemNew(sizeof(NCBISub));
   nsp->ssp = ssp;
   nsp->submittor_key = StringSave(center);
   MemFree(ssp->sub->cit->descr);
   ssp->sub->cit->descr = StringSave(remark);

   sep_list = NULL;
   if (phrap_fname != NULL) {
      sep_list = ReadPhrapFile (fp);
      sep_list = SetPhrapContigOrder (sep_list, contigs);
   } else if (accn_fname != NULL) {
      sep_list = ReadContigList (fp, coordsOnMaster);
   } else {
      while ((sep = FastaToSeqEntry (fp, TRUE)) != NULL) {
         ValNodeLink (&sep_list, sep);
      }
   }

   cumlength = 0;
   index = 0;
   if (accn_fname != NULL) {
    the_entry = AddSeqOnlyToSubmission (
                        nsp,
                        seqname,
                        NULL,
                        accession,
                        0,
                        MOLECULE_CLASS_DNA,
                        MOLECULE_TYPE_GENOMIC,
                        length,
                        TOPOLOGY_LINEAR,
                        STRANDEDNESS_DOUBLE);

      sep = sep_list;
      nextsep = sep->next;
      sep->next = NULL;
      if (the_entry != NULL && the_entry->choice == 1 && sep != NULL && sep->choice == 1) {
         the_bsp = (BioseqPtr)(the_entry->data.ptrvalue);
         bsp = (BioseqPtr)(sep->data.ptrvalue);
         if (the_bsp->repr == Seq_repr_raw) {
            the_bsp->seq_data = BSFree (the_bsp->seq_data);
            the_bsp->repr = bsp->repr;
            the_bsp->seq_data = bsp->seq_data;
            bsp->seq_data = NULL;
            the_bsp->seq_data_type = bsp->seq_data_type;
            the_bsp->seq_ext_type = bsp->seq_ext_type;
            the_bsp->seq_ext = bsp->seq_ext;
            bsp->seq_ext = NULL;
         }
         index++;
      cumlength += bsp->length;
      SeqEntryFree(sep);
      }
   } else if (htgs_phase < 3)
   {
      the_entry = AddDeltaSeqOnlyToSubmission (
                        nsp,
                        seqname,
                        NULL,
                        accession,
                        0,
                        MOLECULE_CLASS_DNA,
                        MOLECULE_TYPE_GENOMIC,
                        length,
                        TOPOLOGY_LINEAR,
                        STRANDEDNESS_DOUBLE);

      sep = sep_list;
      lastwasraw = FALSE;
      while (sep != NULL)
      {
         nextsep = sep->next;
         sep->next = NULL;
         bsp = (BioseqPtr)(sep->data.ptrvalue);
         if (bsp->repr == Seq_repr_raw)
         {
            if (lastwasraw) {
               AddGapToDeltaSeq(nsp, the_entry, 0);
               index++;
            }
            BioseqRawConvert(bsp, Seq_code_iupacna);
            seqbuf = BSMerge((ByteStorePtr)(bsp->seq_data), NULL);
            slp = AddLiteralToDeltaSeq(nsp, the_entry,
               bsp->length);
            AddBasesToLiteral(nsp, slp, seqbuf);
            MemFree(seqbuf);
            lastwasraw = TRUE;
            index++;
         }
         else
         {
            if (bsp->length < 0)
               bsp->length = 0;  /* -1 may be set */
            AddGapToDeltaSeq(nsp, the_entry,
               bsp->length);
            lastwasraw = FALSE;
            index++;
         }
         cumlength += bsp->length;
         RescueSeqGraphs (bsp, index, &rescuedsgps);
         SeqEntryFree(sep);
         sep = nextsep;
      }
   }
   else
   {
    the_entry = AddSeqOnlyToSubmission (
                        nsp,
                        seqname,
                        NULL,
                        accession,
                        0,
                        MOLECULE_CLASS_DNA,
                        MOLECULE_TYPE_GENOMIC,
                        length,
                        TOPOLOGY_LINEAR,
                        STRANDEDNESS_DOUBLE);

      sep = sep_list;
      nextsep = sep->next;
      sep->next = NULL;
      bsp = (BioseqPtr)(sep->data.ptrvalue);
      if (bsp->repr == Seq_repr_raw)
      {
         BioseqRawConvert(bsp, Seq_code_iupacna);
         seqbuf = BSMerge((ByteStorePtr)(bsp->seq_data), NULL);
         AddBasesToBioseq(nsp, the_entry, seqbuf);
         MemFree(seqbuf);
         index++;
      }
      cumlength += bsp->length;
      RescueSeqGraphs (bsp, index, &rescuedsgps);
      SeqEntryFree(sep);
   }

   FileClose(fp);
                             
    /* get data from template: pub, organism, and comment */
   if (IS_Bioseq(oldsep))
   {
      bsp = (BioseqPtr)(oldsep->data.ptrvalue);
      prevpnt = &(bsp->descr);
   }
   else
   {
      bssp = (BioseqSetPtr)(oldsep->data.ptrvalue);
      prevpnt = &(bssp->descr);
   }

   bsp = (BioseqPtr)(the_entry->data.ptrvalue);
   if (bsp != NULL) {
     bsp->length = MAX (cumlength, length);
   }

   for (vnp = *prevpnt; vnp != NULL; vnp = next)
   {
      next = vnp->next;
      if (vnp->choice == Seq_descr_pub
               || (vnp->choice == Seq_descr_comment && temp_comment)
       || ((vnp->choice == Seq_descr_org || vnp->choice == Seq_descr_source)
                && temp_org))
      {
         *prevpnt = next;
         vnp->next = NULL;
         ValNodeLink(&(bsp->descr), vnp);
      }
      else
         prevpnt = &(vnp->next);
   }

   if (comment_fstr != NULL)
   {
      ValNodeCopyStr (&(bsp->descr), Seq_descr_comment, comment_fstr);
      MemFree (comment_fstr);
   } 
   
   SeqEntryFree(oldsep);

   if (organism != NULL)
      AddOrganismToEntryNew(nsp, the_entry, organism, NULL, NULL, NULL,
               NULL, NULL, NULL, NULL);

   AddGenomeToEntry(nsp, the_entry, 1);
   if (clone != NULL)
      AddSubSourceToEntry(nsp, the_entry, 3, clone);
   if (chromosome != NULL)
       AddSubSourceToEntry(nsp, the_entry, 1, chromosome);
   if (clone_lib != NULL)
       AddSubSourceToEntry(nsp, the_entry, 11, clone_lib);
   if (map != NULL)
       AddSubSourceToEntry(nsp, the_entry, 2, map);
   if (strain != NULL)
       AddOrgModToEntry(nsp, the_entry, 2, strain);

   if (title != NULL)
      AddTitleToEntry(nsp, the_entry, title);

   if (htgsDraft) {
      AddGenBankBlockToEntry (nsp, the_entry, NULL, NULL, "HTGS_DRAFT", NULL, NULL);
   }

   if (extra_ac != NULL)
      AddExtraAc2Entry(the_entry, extra_ac);

   AddBiomolToEntry(nsp, the_entry, 1);
   AddTechToEntry(nsp, the_entry, MI_htgs_phase);

   if (bsp != NULL) {
     for (vnp = rescuedsgps; vnp != NULL; vnp = vnp->next) {
       OffsetAndLinkSeqGraph (bsp, (SeqGraphPtr) vnp->data.ptrvalue, (Int2) vnp->choice);
       vnp->data.ptrvalue = NULL;
     }
   }
   rescuedsgps = ValNodeFreeData (rescuedsgps);

   errs = NCBISubValidate(nsp, NULL);

   NCBISubWrite(nsp, myargs[2].strvalue);

   NCBISubFree(nsp);

   return(errs);
}

