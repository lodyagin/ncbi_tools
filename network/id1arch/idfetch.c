/*  idfetch.c
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

* Author Karl Sirotkin
*
$Log: idfetch.c,v $
Revision 1.1  1998/12/28 17:56:29  yaschenk
preparing idfetch to go to production

Revision 1.1  1997/05/29 14:34:07  sirotkin
syncing sampson from mutant for procs. taking source from sampson. this is now current

 * Revision 4.0  1995/07/26  13:55:55  ostell
 * force revision to 4.0
 *
 * Revision 1.3  1995/06/21  14:14:29  kans
 * replaced asn2ff_entrez with SeqEntryToFlat
 *
 * Revision 1.2  1995/05/17  17:59:15  epstein
 * add RCS log revision history
 *
 * Revision 1.1  94/08/11  13:26:31  ostell
 * Initial revision
 * 
 * Revision 1.3  1993/12/02  10:12:41  kans
 * Includes <ncbi.h> instead of <sys/types.h>
 *
 * Revision 1.2  93/11/24  13:25:56  sirotkin
 * First working version
 * 
 * Revision 1.1  93/11/23  16:01:51  sirotkin
 * Initial revision
 * 
   revised by OStell for public use.
 * 
 * Modified by Eugene Yaschenko for ID1 Server
*
*
* RCS Modification History:
* $Log: idfetch.c,v $
* Revision 1.1  1998/12/28 17:56:29  yaschenk
* preparing idfetch to go to production
*
* Revision 1.1  1997/05/29 14:34:07  sirotkin
* syncing sampson from mutant for procs. taking source from sampson. this is now current
*
 * Revision 4.0  1995/07/26  13:55:55  ostell
 * force revision to 4.0
 *
 * Revision 1.3  1995/06/21  14:14:29  kans
 * replaced asn2ff_entrez with SeqEntryToFlat
 *
 * Revision 1.2  1995/05/17  17:59:15  epstein
 * add RCS log revision history
 *
*/
#include <ncbi.h>
#include <objsset.h>
#include <accid1.h>
#include <asn2ff.h>
#include <tofasta.h>
#include <ni_types.h>

Args myargs[] = {
	{"Filename for output ","stdout", NULL,NULL,FALSE,'o',ARG_FILE_OUT, 0.0,0,NULL},
	{"Output type: 1=text asn.1 2=binary asn.1 3=genbank (Seq-entry only) 4=genpept (Seq-entry only) 5=fasta (table for history)",
	  "1", "1", "5", FALSE, 't', ARG_INT, 0.0, 0, NULL } ,
	{"Database to use",NULL,NULL,NULL,TRUE,'d',ARG_STRING,0.0,0,NULL},
	{"Entity number ( retrieval number ) to dump" ,"0","0","99999999",TRUE,'e',ARG_INT,0.0,0,NULL},
        {"Type of lookup:\t\
0 - get Seq-entry\n\t\t\t\
1 - get gi state (output to stderr)\n\t\t\t\
2 - get SeqIds\n\t\t\t\
3 - get gi historyn (sequence change only)\n\t\t\t\
4 - get gi revision history (any change to asn.1)\n", "0","0","4",TRUE,'i',ARG_INT,0.0,0,NULL},
	{"GI id for single Entity to dump" ,"0","1","99999999",TRUE,'g',ARG_INT,0.0,0,NULL},
	{"Maximum complexity of Entity dump (only for -i 0 )" ,"0","0","4",TRUE,'c',ARG_INT,0.0,0,NULL},
 	{"flaTtened SeqId, format: \n		\'type(name,accession,release,version)\'\n			as \'5(HUMHBB)\' or \n		type=accession, or \n		type:number ",
		NULL,NULL,NULL,TRUE,'f',ARG_STRING,0.0,0,NULL},
 	{"Fasta style SeqId ENCLOSED IN QUOTES: lcl|int or str bbs|int bbm|int gb|acc|loc emb|acc|loc pir|acc|name sp|acc|name pat|country|patent|seq gi|int dbj|acc|loc prf|acc|name pdb|entry|chain  ",
	NULL,NULL,NULL,TRUE,'s',ARG_STRING,0.0,0,NULL},
        {"Log file for errors and feedback about IDs assigned", NULL,NULL,NULL,TRUE,'l',ARG_FILE_OUT,
                0.0,0,NULL}
};
int Numarg = sizeof(myargs)/sizeof(myargs[0]);

#define MACRO_SETARG(TAG,P) \
   {\
      P = Nlm_WhichArg (TAG, Numarg, myargs);\
      if ( P < 0){\
         ErrPost(CTX_NCBIIDRETRIEVE,10,\
         "Program error looking for arg %c\n", TAG);\
         has_trouble = TRUE;\
      }\
   }
static Nlm_Int2 Nlm_WhichArg PROTO (( Nlm_Char which, Nlm_Int2 numargs, Nlm_ArgPtr ap));

DataVal Val;

Int2 Main()
{
	AsnIoPtr aip, check_aip = NULL;
	SeqEntryPtr sep = NULL, hold_entry;
	Int2 fileoutarg, logarg, outtypearg,maxplexarg,
		 seqidarg, giarg, fastaarg,infotypearg,entarg,dbarg;
	Boolean has_trouble = FALSE;
	char * msg;
	CharPtr outmode;
	Int4 entity_spec_count = 0;
	Int4 ent = 0;
	Int4 gi = 0;
	Int4 gi_state;
	AsnIoPtr asnout=NULL;
	FILE * fp = NULL;
	Int4 status;
	SeqIdPtr sip,sip_ret;
	ID1SeqHistPtr ishp;
        Char tbuf[40],buf[200];
	

					/* check command line arguments */

	if ( ! GetArgs("idfetch.c",Numarg, myargs))
		return 1;


/********************************************************************
 ****                                                            ****
 ****  Map Args So Can be Accessed in order independent fashion  ****
 ****                                                            ****
 *******************************************************************/

	MACRO_SETARG('o', fileoutarg)
	MACRO_SETARG('t', outtypearg)
        MACRO_SETARG('i', infotypearg)
	MACRO_SETARG('c', maxplexarg)
	MACRO_SETARG('e',entarg)
	MACRO_SETARG('d',dbarg)
	MACRO_SETARG('g',giarg)
	MACRO_SETARG('f',seqidarg)
	MACRO_SETARG('l',logarg)
	MACRO_SETARG('s',fastaarg)

	if (! SeqEntryLoad())
		ErrShow();

	if (myargs[logarg].strvalue != NULL) {
			if (! ErrSetLog (myargs[logarg].strvalue)){
				ErrShow();
				has_trouble = TRUE;
			}else{
				ErrSetOpts (ERR_TEE, ERR_LOG_ON);
			}
	}
	if(myargs[infotypearg].intvalue>1 && (myargs[outtypearg].intvalue == 3 || myargs[outtypearg].intvalue == 4)){
		ErrPostEx(SEV_FATAL,0,0,"-t 3 or -t 4 can be used only with -i 0");
		has_trouble=TRUE;
		goto FATAL;
	}

	if (myargs[giarg].intvalue){
		entity_spec_count ++;
	}
	if (myargs[seqidarg].strvalue){
		entity_spec_count ++;
	}
	if (myargs[fastaarg].strvalue)
		entity_spec_count++;
	
	if (entity_spec_count != 1){
		ErrPostEx(SEV_FATAL,0,0, "One and only one of the -g, -f, -s parameters must be used");
		has_trouble=TRUE;
		goto FATAL;
	}
	if(myargs[infotypearg].intvalue != 1){
		outmode = "w";
		switch (myargs[outtypearg].intvalue)
		{
			case 2:
				outmode = "wb";
			case 1:
				asnout = AsnIoOpen((CharPtr)myargs[fileoutarg].strvalue, outmode);
				if (asnout == NULL)
				{
				 ErrPost(CTX_NCBIIDRETRIEVE,10,
				"Could not open %s for asn output\n",
				    myargs[fileoutarg].strvalue);
				 has_trouble = TRUE;
				}
				break;
			case 3:
			case 4:
			case 5:
				fp = FileOpen((CharPtr)myargs[fileoutarg].strvalue, outmode);
				if (fp == NULL)
				{
				 ErrPost(CTX_NCBIIDRETRIEVE,10,
				"Could not open %s for output\n",
				    myargs[fileoutarg].strvalue);
				 has_trouble = TRUE;
				}
				break;
		}
	}

	if ( has_trouble )
		exit (1);
	NI_SetInterface(eNII_WWWDirect);
	if (!ID1BioseqFetchEnable("idfetch",TRUE)){
		ErrPost(CTX_NCBIIDRETRIEVE,20,
		"Could not open ID1 service");
		exit(1);
	}
	if (myargs[giarg].intvalue){
		gi = myargs[giarg].intvalue;
	}
	else if (myargs[fastaarg].strvalue != NULL)
	{
		sip = SeqIdParse((CharPtr)myargs[fastaarg].strvalue);
		if (sip == NULL)
		{
#ifdef IDFETCH_HTML_OUTPUT
                        fprintf(fp,"<HR><h2>Couldn't parse FASTA format: <I>%s</I></h2>",myargs[fastaarg].strvalue);
			fflush(fp);
#endif
			ErrPostEx(SEV_FATAL,0,0,"Couldn't parse [%s]", myargs[fastaarg].strvalue);
			exit(1);
		}
	}else{
/*  "flaTtened SeqId, format:
            type(name,accession,release,version) or type=accession",
   */
      static CharPtr name = NULL, accession = NULL, release = NULL, version = NULL, number = NULL;
      CharPtr p ;
      int type_int;
      static CharPtr PNTR fields [] = { & name, & accession, & release, & number};      
      Boolean found_equals = FALSE, found_left = FALSE, 
         found_colon = FALSE, flat_seqid_err = FALSE,
         dna_type = FALSE, any_type = FALSE;
      int dex;
      CharPtr sql_where, sql_and, temp_where, temp_and;
		TextSeqIdPtr tsip;

    sip = ValNodeNew(NULL);
		  type_int = atoi(myargs[seqidarg].strvalue);
      for (p = myargs[seqidarg].strvalue; *p; p++ ) {
         if ( *p == '(' || *p == '='  || *p == ':'){  /* ) match */
      
            if ( *p == '('  ){  /* ) match */
               found_left = TRUE;
               if (p == myargs[seqidarg].strvalue){
                  any_type = TRUE;
							ErrPost(CTX_NCBIIDRETRIEVE,10,
							"Sorry, any type is not allowed for ID1service");
							exit(1);
               }else if ( p - myargs[seqidarg].strvalue == 1){
                  if (*myargs[seqidarg].strvalue == '0'){
                     dna_type = TRUE;
							ErrPost(CTX_NCBIIDRETRIEVE,10,
							"Sorry, 0== nucroe3 type is not allowed for ID1service");
							exit(1);
                  }
               }
            }else if ( *p == '=') {
               found_equals = TRUE;
               if (p == myargs[seqidarg].strvalue){
                  any_type = TRUE;
               }else if ( p - myargs[seqidarg].strvalue == 1){
                  if (*myargs[seqidarg].strvalue == '0'){
                     dna_type = TRUE;
                  }
               }
            }else if ( *p == ':'){
               found_colon = TRUE;
               if (p == myargs[seqidarg].strvalue){
                  any_type = TRUE;
               }else if ( p - myargs[seqidarg].strvalue == 1){
                  if (*myargs[seqidarg].strvalue == '0'){
                     dna_type = TRUE;
                  }
               }
            }
            *p = '\0';
            p++;
            break;
         }
      }
      if ( found_left){
         for ( * (fields[0]) = p, dex = 0; *p && dex < 4; p++){
            if ( *p == ',' || *p == ')' ){
                *p = '\0';
               dex ++;
               *(fields[dex]) = p + 1;
            }
         }
      }else if (found_equals){
         accession = p;
      }else if (found_colon){
         number = p;
      }else{
            ErrPost(CTX_NCBIIDRETRIEVE, 10,
            "id1test: could not find \'(\' or \'=\' or \':\' in flattened seqid=%s",myargs[seqidarg].strvalue);  /* ) match */
					exit(1);
      }
		sip -> choice = type_int;
		switch (type_int){
		case SEQID_GIBBSQ :
		case SEQID_GIBBMT :
		case SEQID_GI :
			sip -> data.intvalue = atoi(number);
			break;
		case SEQID_GENBANK : case SEQID_EMBL : case SEQID_DDBJ :
		case SEQID_PIR : case SEQID_SWISSPROT : case SEQID_OTHER :
		case SEQID_PRF :
			tsip = TextSeqIdNew();
			sip -> data.ptrvalue = tsip;
			if (accession)
				if (!*accession)
					accession = NULL;
			if (release)
				if (!*release)
					release = NULL;
			if (name)
				if (!*name)
					name = NULL;
			tsip -> name = StringSave(name);
			tsip -> accession = StringSave(accession);
			tsip -> release = StringSave(release);
			break;
		case SEQID_PATENT : case SEQID_GENERAL : case SEQID_PDB :
		case SEQID_LOCAL :
		ErrPost(CTX_NCBIIDRETRIEVE,30,
		"Sorry, this test program does not support %d patent, general, pdb, or local, try id2asn ", 
			type_int);
			exit(1);
			break;
		}
	}
	if (! gi){
		gi = ID1ArchGIGet (sip);
		if (gi <= 0){
			SeqIdPrint(sip, tbuf, PRINTID_FASTA_SHORT);
#ifdef IDFETCH_HTML_OUTPUT
			fprintf(fp,"<HR><h2>Couldn't find SeqId: <I>%s</I></h2>",tbuf);
			fflush(fp);
#endif
			ErrPostEx(SEV_FATAL, 0,0,"Couldn't find SeqId [%s]", tbuf);
			exit(1);
		}
	}
	switch(myargs[infotypearg].intvalue){
		case 0:
			sep = ID1ArchSeqEntryGet (gi, myargs[dbarg].strvalue, myargs[entarg].intvalue,
					&status,(Int2) myargs[maxplexarg].intvalue);
			if ( !sep){
				switch(status){
				 case 1:
					fprintf(stderr,"Sequence has been withdrawn\n");
					break;
				 case 2:
					fprintf(stderr,"Sequence is not yet available\n");
					break;
				 default:
					fprintf(stderr," Unable to read ASN.1 message\n");
				}
#ifdef IDFETCH_HTML_OUTPUT
				printf("<HR><h2>Sorry, Sequence is not available</h2>");
				fflush(stdout);
#endif
					goto FATAL;
			}
			if (status==3)
				fprintf(stderr," IS DEAD!\n");
			break;
		case 1:
			gi_state = ID1ArcgGIStateGet(gi);
			break;
		case 2:
			sip_ret = ID1ArchSeqIdsGet(gi,asnout);
			break;
		case 3:
			ishp = ID1ArchGIHistGet(gi,FALSE,asnout);
			break;
		case 4:
			ishp = ID1ArchGIHistGet(gi,TRUE,asnout);
			break;
	}
	if(myargs[infotypearg].intvalue == 1){
		Char	buf[200];
		id_print_gi_state(gi_state,buf,sizeof(buf));
		printf("gi= %d, states: %s\n",gi,buf);
	} else {
		switch (myargs[outtypearg].intvalue)
		{
		 case 1:
		 case 2:
			switch(myargs[infotypearg].intvalue){
			 case 0:
				SeqEntryAsnWrite(sep, asnout, NULL);
				break;
			}
			break;
		 case 3:
                        if(!SeqEntryToFlat(sep, fp, GENBANK_FMT, RELEASE_MODE)){
                                ErrPostEx(SEV_WARNING,0,0,
                           "GenBank Format does not exist for this sequence");
				has_trouble=TRUE;
#ifdef IDFETCH_HTML_OUTPUT
                                 fprintf(fp,
                         "<HR><h2>GenBank Format does not exist for this sequence</h2>");
                                fflush(fp);
#endif
			}
			break;
		 case 4:
			if(!SeqEntryToFlat(sep, fp, GENPEPT_FMT, RELEASE_MODE)){
				ErrPostEx(SEV_WARNING,0,0,
			   "GenPept Format does not exist for this sequence");
				has_trouble=TRUE;
#ifdef IDFETCH_HTML_OUTPUT
				 fprintf(fp,
			 "<HR><h2>GenPept Format does not exist for this sequence</h2>");
				fflush(fp);
#endif
			}
			break;
		 case 5:
			switch(myargs[infotypearg].intvalue){
                         case 0:
				SeqEntryToFasta(sep, fp, TRUE);  /* nuc acids */
				SeqEntryToFasta(sep, fp, FALSE); /* proteins */
				break;
			 case 2:
				SeqIdWrite(sip_ret,buf,PRINTID_FASTA_LONG,sizeof(buf) - 1);
				fprintf(fp,"%s\n",buf);
				break;
			 case 3:
			 case 4:
				SeqHistPrintTable(ishp,fp);
				break;
			}
			break;
		}
	}
	SeqEntryFree(sep);
FATAL:
	if(asnout)
		AsnIoClose(asnout);
	if(fp)
		FileClose(fp);
	ID1ArchFini();

	return(has_trouble?1:0);
}

/*****************************************************************************
*
*   Nlm_WhichArg(ap)
*     returns array position for a tag 
*
*****************************************************************************/
static Nlm_Int2 Nlm_WhichArg( Nlm_Char which, Nlm_Int2 numargs, Nlm_ArgPtr ap)
{
   Nlm_Boolean okay = FALSE;
   Nlm_Int2 i;
   Nlm_ArgPtr curarg;
   Nlm_Int2 retval = -1;

   if ((ap == NULL) || (numargs == 0) )
      return okay;

   curarg = ap;                        /* set defaults */

   for (i = 0; i < numargs; i++, curarg++)
   {
      if (curarg->tag == which)
      {
         retval = i;
         break;
      }
   }

   return retval;
}
