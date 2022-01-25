/*
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
* File Name: netblap3.c
*
* Author: Tom Madden
*
* Version Creation Date:   05/8/97
*
* File Description: 
*       Application Programming Interface (API) for BLAST network server
*
* RCS Modification History:
* $Log: netblap3.c,v $
* Revision 1.27  1999/01/12 21:05:58  victorov
* server will now report an error if the ni-queue if full
*
* Revision 1.26  1998/12/09 15:27:05  madden
* Add wordsize
*
* Revision 1.25  1998/11/03 21:46:19  egorov
* Add support for entrez-query and org-name
*
* Revision 1.24  1998/10/26 19:42:57  madden
* Check for NULL matrix or filter_string
*
* Revision 1.23  1998/10/06 18:28:07  egorov
* Fix MT problem with dispatcher_lock
*
* Revision 1.22  1998/09/22 16:12:53  egorov
* Use BlastErrorPrintExtra instead of BlastErrorPrint to redirect error messages both to log and to program output file
*
* Revision 1.21  1998/09/01 20:17:03  madden
*  Fixed uninitialzed problem in BlastNetBioseqFetchDisable, changed prototype
*
* Revision 1.20  1998/08/14 17:43:26  madden
* non-NULL f_order and g_order in formatting
*
* Revision 1.19  1998/08/14 16:04:48  egorov
* Create TLS for BlastNetFetchStruct to make it multi-thread safe.  Using some of global variables is commented out
*
* Revision 1.18  1998/07/28 16:57:51  egorov
* Make StringSave() for some CharPtr in assignments to save memory and avoid memory crashes
*
* Revision 1.16  1998/06/08 21:58:36  madden
* Check for NULL Bioseq in FetchFunc, print out Error message if NULL
*
* Revision 1.15  1998/05/08 21:40:15  vakatov
* fixed a tiny type cast(one more *)
*
* Revision 1.14  1998/05/08 20:56:22  madden
* Fix PC compile warnings, rename callback
*
* Revision 1.13  1998/05/08 01:08:09  madden
* Removed unused variables
*
* Revision 1.12  1998/04/24 18:35:50  madden
* Do not send Seq-descr to server
*
* Revision 1.11  1998/04/24 16:01:41  egorov
* Remove BlastErrorPrint from BlastSeqLocNetCore
*
* Revision 1.10  1998/04/23 14:19:31  egorov
* Allow other_returns in BlastSeqLocNet to be NULL
*
* Revision 1.9  1998/04/22 18:10:06  egorov
* Add support for SeqLoc to blastcl3
*
* Revision 1.8  1998/04/17 19:24:00  madden
* Transmit expect value, hitlist_size, and genetic codes to server
*
* Revision 1.7  1998/04/16 19:35:31  madden
* Added Int4Ptr arg to TraditionalBlastReport specifying the numbers of hits
*
* Revision 1.6  1998/04/09 16:16:33  madden
* Added BlastNetFetchCompare function
*
* Revision 1.5  1997/12/16 19:10:51  madden
* Fixed Codecenter errors
*
* Revision 1.4  1997/12/01 22:05:51  madden
* Changed call to BlastValidateOptionsEx
*
* Revision 1.3  1997/11/28 18:17:39  madden
* Changes for multiple db searches
*
* Revision 1.1  1997/10/08 19:27:20  madden
* Network support for gapped blast
*
*/

#include <ncbinet.h>
#include <ncbithr.h>
#include <txalign.h>
#include <jzmisc.h>
#include <blastpat.h>
#include <netblap3.h>


#define BLAST_SERVER_RETRIES  2

#define BLASTNET_SHORT_BUFLEN 25
#define BLASTNET_BUF_SIZE 255
#define BLASTNET_INIT 0
#define BLASTNET_READY 1
#define BLASTNET_DISABLE 2

typedef struct _blastnetbioseqfetch {
	Uint1 BlastNetFetchState;
	CharPtr dbname;	/* Name of the database. */
	BlastNet3Hptr bl3hp;	/* Network connection. */
	Uint2 ctr;
	Boolean is_prot; /* Is it a protein or not. */
} BlastNetFetchStruct, PNTR BlastNetFetchStructPtr;

static TNlmTls blastnetfetch_tls = NULL;

static Int2 num_attached = 0;
static NI_DispatcherPtr dispatcher = NULL;
/* Mutex for dispatcher. */
TNlmMutex dispatcher_lock = NULL;

static TNlmMutex formating_mutex; /* Mutex to regulate formating in TraditionalBlastReport */

static Boolean BlastInitEx PROTO((CharPtr program_name, BlastNet3Hptr bl3hp, BlastResponsePtr PNTR respp));
 
static Boolean SubmitRequest PROTO((BlastNet3Hptr bl3hptr, BlastRequestPtr blreqp, BlastResponsePtr PNTR response, NetProgressCallback callback));
static Boolean RealSubmitRequest PROTO((BlastNet3Hptr bl3hptr, BlastRequestPtr blreqp, BlastResponsePtr PNTR response, NetProgressCallback callback));

/* error_occurred and old_error_hook should not be accessed directly.
The functions BlastSetErrorStatus, BlastGetErrorStatus, and BlastSetErrorHook
will change these variables. */
static Boolean error_occurred;
static ErrHookProc old_error_hook;

/*
	The next five functions set an error hook and detect
	whether an error occurred (i.e., contact with the
	server was lost).

	BlastSetErrorHook should be called first, the function
	BlastErrHookProc is set as the "handler"; BlastSetErrorStatus
	should then be called to set the error status to FALSE; 
	BlastGetErrorStatus should be called to determine if an error 
	occurred; and then BlastResetOldHook should be called to restore 
	the original hook, in case BLAST is called from within another
	application that uses this (original) hook.

	BlastSetErrorStatus is also called, to set the error status
	sometimes if a BlastAsnRead fails.  This is taken as evidence
	of an error, even if none is reported.
*/

static int LIBCALLBACK
BlastErrHookProc(const ErrDesc *err)

{
	error_occurred = TRUE;
	return 1;
}

static void
BlastSetErrorHook(void)

{
	old_error_hook = Nlm_ErrSetHandler(BlastErrHookProc);
	return;
}

static void
BlastSetErrorStatus(Boolean status)

{
	error_occurred = status;
}

static Boolean
BlastGetErrorStatus(void)

{
	return error_occurred;
}

static void
BlastResetOldHook(void)

{
	Nlm_ErrSetHandler(old_error_hook);
	return;
}

/*
	The following functions fill a the Error user string with
	text to identify BLAST and the entry being worked on.
	The SeqIdPtr is used to make a FASTA id, which is appended
	to string.

	A Uint1 is returned, which allows Nlm_ErrUserDelete to delete
	this error string when it's done.
*/

static Uint1
BlastSetUserErrorString(CharPtr string, SeqIdPtr sip)

{
	Char buffer[2*BLASTNET_SHORT_BUFLEN], textid[BLASTNET_SHORT_BUFLEN];
	CharPtr buf_start, ptr;
	Int2 length;

	ptr = buf_start = &buffer[0];

	StringNCpy(ptr, string, BLASTNET_SHORT_BUFLEN);
	if (sip != NULL)
	{
	    length = StringLen(string);
	    if (length > BLASTNET_SHORT_BUFLEN)
		length = BLASTNET_SHORT_BUFLEN;

	    ptr += length;

    	    SeqIdWrite(sip, textid, PRINTID_FASTA_LONG, BLASTNET_SHORT_BUFLEN-1);
	    StringNCpy(ptr, textid, BLASTNET_SHORT_BUFLEN-1);
	}
	return Nlm_ErrUserInstall (buf_start, 0);
}

static void
BlastDeleteUserErrorString(Uint1 err_id)

{
	Nlm_ErrUserDelete(err_id);
	return;
}

/*****************************************************************************
*
*   NetInit ()
*
*****************************************************************************/

static Boolean NetInit(void)
{
	Boolean retval = FALSE;

	NlmMutexLockEx(&dispatcher_lock);

	if(num_attached++ > 0)
	{
	    NlmMutexUnlock(dispatcher_lock);
	    return TRUE;
	}
        
    	dispatcher = NI_GenericInit(NULL, NULL, TRUE, NULL, 0);

	if (dispatcher)
		retval = TRUE;

	NlmMutexUnlock(dispatcher_lock);

	return retval;
}


/*****************************************************************************
*
*   ForceNetInit ()
*
*****************************************************************************/

static Boolean ForceNetInit(void)
{
	Boolean retval;

	NlmMutexLockEx(&dispatcher_lock);

	num_attached = 0; /* force re-attempt to contact dispatcher */
	retval = NetInit();

	NlmMutexUnlock(dispatcher_lock);

	return retval;
}

/*****************************************************************************
*
*   NetFini ()

	BlastNet3Hptr bl3hp: was returned by BlastInit
	Boolean deallocate: should BlastNet3Hptr be deallocated (or will it be reused for
		a reconnection).
*
*****************************************************************************/

static Boolean NetFini(BlastNet3Hptr bl3hp, Boolean deallocate)
{
    NlmMutexLockEx(&dispatcher_lock);

	if (num_attached > 0)
		num_attached--;

	if (bl3hp != NULL)
	{
		NI_ServiceDisconnect(bl3hp->svcp);
		if (deallocate)
			bl3hp = MemFree(bl3hp);

		if (num_attached == 0)
		{	/* Disconnect if last service to dispatcher. */
			NI_EndServices (dispatcher);
			dispatcher = NULL;
		}
	}

	NlmMutexUnlock(dispatcher_lock);

    return TRUE;
}


/*****************************************************************************
*
*   GenericReestablishNet ()
*
*****************************************************************************/

static Boolean GenericReestablishNet(CharPtr svcName, Boolean showErrs, BlastNet3Hptr bl3hp)
{
    Handle mon = NULL;
    Boolean retval;
    CharPtr buf;

    buf = MemNew(2 * StrLen(svcName) + 60);

    if (showErrs) {
        sprintf (buf, "Re-establishing %s Service", svcName);
        mon = MonitorStrNew(buf, 40);
        sprintf (buf, "Requesting %s service", svcName);
        MonitorStrValue(mon, buf);
    }

    NetFini(bl3hp, FALSE);

    retval = BlastInitEx(NULL, bl3hp, NULL);

    if (!retval)
    {
        sprintf (buf, "%s get failed; re-contacting dispatcher", svcName);
        MonitorStrValue(mon, buf);
        retval = FALSE;
        if (ForceNetInit())
        { /* successfully established contact w/dispatcher */
            sprintf (buf, "%s get failed; re-requesting %s service",
                     svcName, svcName);
            MonitorStrValue(mon, buf);
	    retval = BlastInitEx(NULL, bl3hp, NULL);
        }
        else {
            ErrPost(CTX_UNKNOWN, 1, "Unable to re-contact dispatcher");
            if (showErrs) {
                ErrShow();
            }
        }
    }

    MonitorFree(mon);

    if (! retval )
    {
        sprintf (buf, "Unable to re-establish %s service", svcName);
        ErrPost(CTX_UNKNOWN, 1, buf);
        if (showErrs) {
            ErrShow();
        }
    }

    MemFree(buf);
    return retval;
}

/*****************************************************************************
*
*   ReestablishNetBlast ()
*
*****************************************************************************/

static Boolean ReestablishNetBlast(BlastNet3Hptr bl3hp)
{
    return GenericReestablishNet("NETBLAST", TRUE, bl3hp);
}

/*****************************************************************************
*
*   BlastInit ()
*
*****************************************************************************/

static Boolean BlastInitEx (CharPtr program_name, BlastNet3Hptr bl3hp, BlastResponsePtr PNTR respp)

{
	BlastRequestPtr request;
	BlastResponsePtr response;
	NI_HandPtr svcp = NULL;

	if (bl3hp == NULL)
	{
		return FALSE;
	}


	if (!NetInit())
    	{
       	 	ErrPostEx(SEV_ERROR, 0, 0, "NI_ServiceGet [%s] (%s)", ni_errlist[ni_errno], ni_errtext);
		return FALSE;
	}
	
    	svcp = NI_GenericGetService(dispatcher, NULL, "NETBLAST", "blast3", FALSE);
    	if (svcp == NULL)
    	{
       	 	ErrPostEx(SEV_ERROR, 0, 0, "NI_ServiceGet [%s] (%s)", ni_errlist[ni_errno], ni_errtext);
       	 	BlastFini(NULL);
       	 	return FALSE;
    	}

	bl3hp->svcp = svcp;


	request = ValNodeNew(NULL);
	request->choice = BlastRequest_init;
	SubmitRequest(bl3hp, request, &response, NULL);
    	BlastRequestFree (request);
	if (respp)
		*respp = response;
	else
    		BlastResponseFree (response);

	return TRUE;

}

/*****************************************************************************
*
*   BlastInit ()
*
*****************************************************************************/

Boolean LIBCALL
BlastInit (CharPtr program_name, BlastNet3Hptr PNTR bl3hpp, BlastResponsePtr PNTR respp)

{
	BlastNet3Hptr bl3hp_pri;

	if (bl3hpp == NULL)
	{
       	 	ErrPostEx(SEV_ERROR, 0, 0, "BlastInitMT, BlastNet3Hptr PNTR is NULL");
		return FALSE;
	}

	bl3hp_pri = MemNew(sizeof(BlastNet3Hptr));
	if (bl3hp_pri == NULL)
		return FALSE;

	*bl3hpp = bl3hp_pri;

	return BlastInitEx(program_name, bl3hp_pri, respp);
}

/*****************************************************************************
*
*   BlastFini ()
*
*****************************************************************************/

static Boolean s_BlastFini (BlastNet3Hptr bl3hptr)

{
	Boolean retval = TRUE;
	BlastRequestPtr request;
	BlastResponsePtr response;


	if (bl3hptr == NULL || bl3hptr->svcp == NULL)
		return FALSE;

        request = ValNodeNew(NULL);
        request->choice = BlastRequest_fini;
	SubmitRequest(bl3hptr, request, &response, NULL);
        BlastRequestFree (request);
        BlastResponseFree (response);

	NetFini(bl3hptr, TRUE);
	return retval;
}

/* the only thing done here is to suppress errors */

Boolean LIBCALL 
BlastFini (BlastNet3Hptr bl3hptr)

{
    short erract;
    ErrDesc err;
    Boolean retval;

    ErrGetOpts(&erract, NULL);
    ErrSetOpts(ERR_IGNORE, 0);
    ErrFetch(&err);

    retval = s_BlastFini(bl3hptr);

    ErrSetOpts(erract, 0);
    ErrFetch(&err);

    return retval;
}

BlastParametersPtr LIBCALL
BlastOptionsToParameters (BLAST_OptionsBlkPtr options)

{
	BlastParametersPtr parameters;


	if (options == NULL)
	{
		return NULL;
	}

	parameters = BlastParametersNew();
	parameters->gapped_alignment = options->gapped_calculation;
	if (options->cutoff_s == 0)
	{
		parameters->Cutoff_cutoff = ValNodeAddFloat(NULL, Cutoff_cutoff_evalue, options->expect_value);
	}
	else
	{
		parameters->Cutoff_cutoff = ValNodeAddInt(NULL, Cutoff_cutoff_score, options->cutoff_s);
	}
	if (options->cutoff_s2 == 0)
	{
		parameters->Cutoff2_cutoff2 = ValNodeAddFloat(NULL, Cutoff2_cutoff2_evalue, options->e2);
	}
	else
	{
		parameters->Cutoff2_cutoff2 = ValNodeAddInt(NULL, Cutoff2_cutoff2_score, options->cutoff_s2);
	}
	parameters->hitlist_size = options->hitlist_size;
	parameters->first_threshold = options->threshold_first;
	parameters->second_threshold = options->threshold_second;
	parameters->nucl_penalty = options->penalty;
	parameters->nucl_reward = options->reward;
	parameters->gap_open = options->gap_open;
	parameters->gap_extend = options->gap_extend;
	parameters->genetic_code = options->genetic_code;
	parameters->db_genetic_code = options->db_genetic_code;
	parameters->low_complexity_filtering = options->filter;
	parameters->ethresh = options->ethresh;
        parameters->max_num_passes = options->maxNumPasses;
        parameters->pseudo_count_const = options->pseudoCountConst;

        parameters->gifile = StringSave(options->gifile);
        parameters->gilist = options->gilist;
        parameters->matrix = StringSave(options->matrix);
        parameters->filter_string = StringSave(options->filter_string);
	parameters->entrez_query = StringSave(options->entrez_query);
        parameters->word_size = options->wordsize;
        
	return parameters;
}

/*
	Translates the BlastDbinfoPtr into TxDfDbInfoPtr.
*/

TxDfDbInfoPtr LIBCALL
NetDbinfo2TxDbinfo(BlastDbinfoPtr net_dbinfo)


{
	TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;

	dbinfo_head = NULL;
	while (net_dbinfo)
	{
		dbinfo = TxDfDbInfoNew(dbinfo);
		dbinfo->is_protein = net_dbinfo->is_protein;
		dbinfo->name = StringSave(net_dbinfo->name);
		dbinfo->definition = StringSave(net_dbinfo->definition);
		dbinfo->date = StringSave(net_dbinfo->date);
		dbinfo->total_length = net_dbinfo->total_length;
		dbinfo->number_seqs = net_dbinfo->number_seqs;
		if (dbinfo_head == NULL)
			dbinfo_head = dbinfo;
		net_dbinfo = net_dbinfo->next;
	}

	return dbinfo_head;
}

BlastNet3BlockPtr LIBCALL
BlastNet3BlockDestruct(BlastNet3BlockPtr blnet)

{
	BlastResponsePtr response;
	if (blnet == NULL)
		return NULL;

	BlastParametersFree(blnet->parameters);
	/* individual elements freed elsewhere. */
	response = blnet->response;
	while(response)
	{
		response->data.ptrvalue = NULL;
		response = response->next;
	}
	BlastResponseFree(blnet->response);
	MemFree(blnet->dbname);

	return NULL;
}

/*
	Creates a new BlastNet3BlockNew, used for searching.
*/
BlastNet3BlockPtr LIBCALL
BlastNet3BlockNew(CharPtr program, CharPtr dbname)

{
	BlastNet3BlockPtr retval=NULL;

	retval = MemNew(sizeof(BlastNet3Block));

	if (retval)
	{
		if (StringICmp("blastn", program) == 0)
			retval->prog_type = Blast_search_program_blastn;
		else if (StringICmp("blastp", program) == 0)
			retval->prog_type = Blast_search_program_blastp;
		else if (StringICmp("blastx", program) == 0)
			retval->prog_type = Blast_search_program_blastx;
		else if (StringICmp("tblastn", program) == 0)
			retval->prog_type = Blast_search_program_tblastn;
		else if (StringICmp("tblastx", program) == 0)
			retval->prog_type = Blast_search_program_tblastx;
		if (dbname)
			retval->dbname = StringSave(dbname);
	}

	return retval;
}

static Boolean
QueryIsProteinFromType(Uint2 type)
{
	switch (type)
	{
		case Blast_search_program_blastn:
			return FALSE;
		case Blast_search_program_blastp:
			return TRUE;
		case Blast_search_program_blastx:
			return FALSE;
		case Blast_search_program_tblastn:
			return TRUE;
		case Blast_search_program_tblastx:
			return FALSE;
		default: 
			return FALSE;
	}
}

static BlastResponsePtr GetResponsePtr(BlastResponsePtr response, Nlm_Uint1 choice)

{
	while (response)
	{
		if (response->choice == choice)
		{
			break;
		}
		response = response->next;
	}

	return response;
}


BlastDbinfoPtr LIBCALL
BlastRequestDbInfo (BlastNet3Hptr bl3hp, CharPtr database, Boolean is_prot)

{
	BlastRequestPtr request=NULL;
	BlastResponsePtr response=NULL;
	BlastDbinfoPtr dbinfo=NULL;
	BlastDbinfoGetPtr dbinfo_get;


	dbinfo_get = BlastDbinfoGetNew();
	dbinfo_get->name = database;

	if (is_prot)
		dbinfo_get->type = Blast_dbinfo_get_type_protein;
	else
		dbinfo_get->type = Blast_dbinfo_get_type_nucleotide;

	ValNodeAddPointer(&request, BlastRequest_db_info_specific, dbinfo_get);
	SubmitRequest(bl3hp, request, &response, NULL);
	response = GetResponsePtr(response, BlastResponse_db_info_specific);

	if (response)
		dbinfo = response->data.ptrvalue;

	dbinfo_get->name = NULL;	/* Not owned by this function. */
        BlastRequestFree(request);

	return dbinfo;
}

BlastDbinfoPtr LIBCALL
BlastGetDbInfo (BlastNet3BlockPtr blnet3blkptr)

{
	BlastResponsePtr response;
	BlastDbinfoPtr dbinfo=NULL, dbinfo_head, last;

	last = NULL;
	dbinfo = NULL;
	dbinfo_head = NULL;
	response = blnet3blkptr->response;
	while (response)
	{
		response = GetResponsePtr(response, BlastResponse_db_info_specific);
		if (response)
		{
			last = dbinfo;
			dbinfo = response->data.ptrvalue;
			if (dbinfo_head == NULL)
				dbinfo_head = dbinfo;
			if (last)
				last->next = dbinfo;
			response = response->next;	
		}
	}

	return dbinfo_head;
}

BlastMatrixPtr LIBCALL
NetBlastGetMatrix(BlastNet3BlockPtr blnet3blkptr)

{
	BlastResponsePtr response;
	BlastMatrixPtr matrix=NULL;
	

	response = GetResponsePtr(blnet3blkptr->response, BlastResponse_matrix);

	if (response)
		matrix = response->data.ptrvalue;

	return matrix;
}

CharPtr LIBCALL
BlastGetParameterBuffer (BlastNet3BlockPtr blnet3blkptr)

{
	BlastResponsePtr response;
	CharPtr buffer=NULL;

	response = GetResponsePtr(blnet3blkptr->response, BlastResponse_parameters);

	if (response)
		buffer = response->data.ptrvalue;

	return buffer;
}

/*
	Find the BlastKABlkPtr of the type (gapped or ungapped) specified.
*/
BlastKABlkPtr LIBCALL
BlastGetKaParams (BlastNet3BlockPtr blnet3blkptr, Boolean gapped)

{
	BlastResponsePtr response;
	BlastKABlkPtr kablk=NULL;


	response = blnet3blkptr->response;
	while (response)
	{
		response = GetResponsePtr(response, BlastResponse_kablk);
		if (response)
		{
			kablk = response->data.ptrvalue;
			if (kablk->gapped == gapped)
				break;
			response = response->next;	
		}
	}

	return kablk;
}

/*
	Converts 'standard' BLAST matrix to network matrix. 
*/
BlastMatrixPtr LIBCALL
BlastMatrixToBlastNetMatrix(BLAST_MatrixPtr matrix)

{
	BlastMatrixPtr net_matrix;
	Int4 index1, index2;
	ValNodePtr vnp=NULL;

	if (matrix == NULL)
		return NULL;

	net_matrix = MemNew(sizeof(BlastMatrix));

	net_matrix->is_protein = matrix->is_prot;
	net_matrix->name = StringSave(matrix->name);
	net_matrix->row_length = matrix->rows;
        net_matrix->column_length = matrix->columns;
	net_matrix->karlinK = matrix->karlinK;

	for (index1=0; index1<matrix->rows; index1++)
	{
		for (index2=0; index2<matrix->columns; index2++)
		{
			ValNodeAddInt(&vnp, 0, matrix->matrix[index1][index2]);
		}
	}
	net_matrix->scores = vnp;

	return net_matrix;

}

/*
	Coverts the 'network' matrix to a 'tools' matrix.
*/
BLAST_MatrixPtr LIBCALL
BlastNetMatrixToBlastMatrix (BlastMatrixPtr net_matrix)

{
	BLAST_MatrixPtr blast_matrix;
	Int4 index1, index2;
	Int4Ptr PNTR matrix;
	ValNodePtr vnp;

	if (net_matrix == NULL)
		return NULL;

	blast_matrix = MemNew(sizeof(BLAST_Matrix));

	blast_matrix->is_prot = net_matrix->is_protein;
	blast_matrix->name = StringSave(net_matrix->name);
	blast_matrix->rows = net_matrix->row_length;
        blast_matrix->columns = net_matrix->column_length;
	blast_matrix->karlinK = net_matrix->karlinK;

	vnp = net_matrix->scores;
	matrix = (Int4Ptr PNTR) MemNew(blast_matrix->rows*sizeof(Int4Ptr));
	for (index1=0; index1<blast_matrix->rows; index1++)
	{
		matrix[index1] = MemNew(blast_matrix->columns*sizeof(Int4));
		for (index2=0; index2<blast_matrix->columns; index2++)
		{
			matrix[index1][index2] = (Int4) vnp->data.intvalue;
			vnp = vnp->next;
		}
	}
	blast_matrix->matrix = matrix;

	return blast_matrix;
}

ValNodePtr LIBCALL
BlastGetMaskedLoc (BlastNet3BlockPtr blnet3blkptr)

{
	BlastResponsePtr response;
	BlastMaskPtr blast_mask;
	ValNodePtr mask_loc=NULL;

	response = blnet3blkptr->response;

	while (response)
	{
		response = GetResponsePtr(response, BlastResponse_mask);
		if (response)
		{
			blast_mask = response->data.ptrvalue;
			ValNodeAddPointer(&mask_loc, blast_mask->frame, blast_mask->location);
			response = response->next;	
		}
	}

	return mask_loc;
}

static BioseqPtr 
PrivateBlastGetBioseq(BlastNet3Hptr bl3hptr, CharPtr database, SeqIdPtr sip, Boolean is_prot)

{
	BioseqPtr bsp=NULL;
	BlastSeqIdPtr blast_sip;
	BlastRequestPtr request = NULL;
	BlastResponsePtr response;

	blast_sip = BlastSeqIdNew();
	blast_sip->is_protein = is_prot;
	blast_sip->database = database;
	blast_sip->id = sip;
	ValNodeAddPointer(&request, BlastRequest_db_seq_get, blast_sip);
	SubmitRequest(bl3hptr, request, &response, NULL);
	
	response = GetResponsePtr(response, BlastResponse_db_seq_get);
	if (response)
	{
		bsp = response->data.ptrvalue;
		response->data.ptrvalue = NULL;
		BlastResponseFree(response);
	}

	/* These two were not allocated here. */
	blast_sip->id = NULL;
	blast_sip->database = NULL;
	BlastRequestFree(request);

	return bsp;
}

BioseqPtr LIBCALL
BlastGetBioseq(BlastNet3BlockPtr blnet3blkptr, SeqIdPtr sip)

{
	BioseqPtr bsp;
	
	bsp = PrivateBlastGetBioseq(blnet3blkptr->bl3hptr, blnet3blkptr->dbname, sip, QueryIsProteinFromType(blnet3blkptr->prog_type));

	return bsp;
}

/*
	Basic functions to submit BLAST runs.

*/

SeqAlignPtr LIBCALL
BlastBioseq (BlastNet3BlockPtr blnet3blkptr, ValNodePtr *error_returns)

{
	BlastRequestPtr request = NULL;
	BlastSearchPtr search = NULL;
	BlastResponsePtr response;
	ValNodePtr node;
	SeqAlignPtr seqalign=NULL;
	Uint1 err_id;

#if 0        
	err_id = BlastSetUserErrorString("netblast:", blnet3blkptr->bsp->id);
#endif        

	search = BlastSearchNew();
	search->program = blnet3blkptr->prog_type;
/*
	search->query = (struct struct_Bioseq *) blnet3blkptr->bsp;
*/
	search->query = blnet3blkptr->bsp;
	search->database = blnet3blkptr->dbname;
	search->parameters = blnet3blkptr->parameters;
	search->mask = blnet3blkptr->mask;
	search->matrix = BlastMatrixToBlastNetMatrix(blnet3blkptr->blast_matrix);
	ValNodeAddPointer(&request, BlastRequest_search, search);
	SubmitRequest(blnet3blkptr->bl3hptr, request, &response, blnet3blkptr->callback);
	
	blnet3blkptr->response = response;
	node = GetResponsePtr(response, BlastResponse_alignment);
	if (node)
		seqalign = node->data.ptrvalue;
	*error_returns = GetResponsePtr(response, BlastResponse_error);

	/* These four are not allocated here. */
	search->query = NULL;
	search->database = NULL;
	search->parameters = NULL;
	search->mask = NULL;
        BlastRequestFree(request);
#if 0        
	BlastDeleteUserErrorString(err_id);
#endif        

	return seqalign;
}

CharPtr LIBCALL 
Blast3GetMotd(BlastNet3Hptr bl3hptr)

{
	BlastResponsePtr response=NULL;
	BlastRequestPtr request=NULL;
	CharPtr string;
	
        request = ValNodeNew(NULL);
        request->choice = BlastRequest_motd;
        SubmitRequest(bl3hptr, request, &response, NULL);
        BlastRequestFree (request);
	if (response == NULL || response->choice != BlastResponse_motd)
		return NULL;

	string = StringSave(response->data.ptrvalue);

	BlastResponseFree(response);

	return string;
}

static Boolean
SubmitRequest(BlastNet3Hptr bl3hptr, BlastRequestPtr blreqp, BlastResponsePtr PNTR response, NetProgressCallback callback)

{
    Boolean status;
    Int2 index;

#if 0    
    BlastSetErrorHook();
#endif    
    for(index=0; index<BLAST_SERVER_RETRIES; index++)
    {
#if 0    
	BlastSetErrorStatus(FALSE);
#endif    

	status = RealSubmitRequest(bl3hptr, blreqp, response, callback);

#if 0    
	if(status == TRUE && BlastGetErrorStatus() == FALSE)
		break;
#else
        if (status == TRUE)
            break;
#endif    

 	ReestablishNetBlast(bl3hptr);
    }
#if 0    
    BlastResetOldHook();
#endif    

    return status;

}

/*
	Status returned indicates a (potential) error.  if not done, this is an error.
*/
#define PRINT_DEBUG_ASN 0
#if PRINT_DEBUG_ASN	
static	tmpi = 0;
#endif 

static Boolean
RealSubmitRequest(BlastNet3Hptr bl3hptr, BlastRequestPtr blreqp, BlastResponsePtr PNTR response, NetProgressCallback callback)

{
	AsnIoPtr asnin, asnout;
	Boolean cancel, done;
	BlastResponsePtr bllist, blresp, head;
#if PRINT_DEBUG_ASN	
        AsnIoPtr asnout_test;
	Char buf[1024];
#endif

	if (bl3hptr == NULL)
		return FALSE;

	if (response == NULL)
		return FALSE;

	asnout = bl3hptr->svcp->waip;
	asnin = bl3hptr->svcp->raip;

#if PRINT_DEBUG_ASN	
        sprintf(buf, "request.%d.out", ++tmpi);
        asnout_test = AsnIoOpen(buf, "w");
	BlastRequestAsnWrite(blreqp, asnout_test, NULL);
#endif	

	done = BlastRequestAsnWrite(blreqp, asnout, NULL);
	if (done == FALSE)
		return FALSE;

#if PRINT_DEBUG_ASN	
        AsnIoReset(asnout_test);
        AsnIoClose(asnout_test);   
#endif	
	AsnIoReset(asnout);

	head = NULL;
	done = FALSE;
	while (!done && (blresp = BlastResponseAsnRead(asnin, NULL)) != NULL)
	{
		switch (blresp->choice)
		{
		     case BlastResponse_done:
			done = TRUE;
			blresp = BlastResponseFree(blresp);
			break;
		    
		     case BlastResponse_queued:
		     case BlastResponse_start:
		     case BlastResponse_progress:
			if (callback)
			{
				if (callback(blresp, &cancel) == TRUE)
				{
					blresp = BlastResponseFree(blresp);
				}
				else
				{
					blresp = BlastResponseFree(blresp);
					done = TRUE;
				}
			}
			break;

		     case BlastResponse_init:
		     case BlastResponse_motd:
		     case BlastResponse_fini:
			done = TRUE;

		     default:
			if (head == NULL)
			{
				head = blresp;
				bllist = blresp;
			}
			else
			{
				bllist->next = blresp;
				bllist = bllist->next;
			}
			break;
		}
	}
	*response = head;
	return done;
}


static void BlastNetFetchCleanup (TNlmTls tls, VoidPtr ptr)
{
	BlastNetFetchStructPtr bnfsp = (BlastNetFetchStructPtr) ptr;

	MemFree(bnfsp);
	return;
}

/*
	Checks the chain of BlastNetFetchStructPtr's for one
	which belongs to the calling thread. If none is found,
	NULL isreturned; otherwise the BlastNetFetchStructPtr is
	returned.
*/
static BlastNetFetchStructPtr
BlastNetFindFetchStruct(BlastNet3Hptr bl3hp, CharPtr dbname, Boolean is_na)

{
	BlastNetFetchStructPtr	bnfsp = NULL;

	if (NlmTlsGetValue(blastnetfetch_tls, (VoidPtr *)(&bnfsp)))
	{
		if (bnfsp == NULL)
		{
			bnfsp = MemNew(sizeof(BlastNetFetchStruct));
			bnfsp->dbname = StringSave(dbname);
			bnfsp->is_prot = (is_na == TRUE) ? FALSE : TRUE;
			bnfsp->bl3hp = bl3hp;
			bnfsp->BlastNetFetchState = BLASTNET_INIT;
			NlmTlsSetValue(&blastnetfetch_tls, bnfsp, 
				BlastNetFetchCleanup);
		}
	}

	return bnfsp;
}

static Boolean BlastNetInit(BlastNetFetchStructPtr bnfsp)

{
	return TRUE;

}


/**********************************************************************

	Fetches the Bioseq, based on the ordinal number of the
	sequence in the database.

************************************************************************/

static Int2 LIBCALLBACK BlastNetBioseqFetchFunc(Pointer data)
{
	Boolean status;
	Char buffer[64];
	OMProcControlPtr ompcp;
        ObjMgrProcPtr ompp;
	OMUserDataPtr omdp;
	SeqIdPtr sip, best_id;
	BlastNetFetchStructPtr blfsp;
	SeqEntryPtr sep;
	BioseqPtr bsp, core_bsp;

	ompcp = (OMProcControlPtr)data;
        ompp = ompcp->proc;

	blfsp = BlastNetFindFetchStruct(NULL, NULL, FALSE);

	if (blfsp->BlastNetFetchState == BLASTNET_INIT)
	{
		status = BlastNetInit(blfsp);
		if (status == FALSE)
			return OM_MSG_RET_OK;
		blfsp->BlastNetFetchState = BLASTNET_READY;
	}

	sip = (SeqIdPtr) (ompcp->input_data);

	best_id = SeqIdFindBest(sip, SEQID_GI);

	if (best_id == NULL)
	{
		core_bsp = BioseqFindCore(sip);
		if (core_bsp)
			best_id = SeqIdFindBest(core_bsp->id, SEQID_GI);
	}

	if (best_id == NULL)
		return OM_MSG_RET_OK;

        
	/* A BioseqPtr is returned by this function. */
	bsp = PrivateBlastGetBioseq(blfsp->bl3hp, blfsp->dbname, best_id, blfsp->is_prot);
        if (bsp == NULL)
        {
                SeqIdWrite(best_id, buffer, PRINTID_FASTA_LONG, sizeof(buffer));
                ErrPost(CTX_UNKNOWN, 1, "Unable to retrieve %s", buffer);
                return OM_MSG_RET_ERROR;
        }
	sep = SeqEntryNew();
	sep->choice = 1;
	sep->data.ptrvalue = bsp;
	SeqMgrSeqEntry(SM_BIOSEQ, (Pointer)bsp, sep);
	ompcp->output_data = (Pointer)bsp;
	ompcp->output_entityID = ObjMgrGetEntityIDForChoice(sep);
	omdp = ObjMgrAddUserData(ompcp->output_entityID, ompp->procid, OMPROC_FETCH, 0);

	return OM_MSG_RET_DONE;
}

/*
	Compare BlastNetFetchStructPtr structures with some traits.  TRUE if identical, otherwise
	FALSE.
*/
static Boolean
BlastNetFetchCompare (BlastNetFetchStructPtr blfsp1, BlastNet3Hptr bl3hp, CharPtr dbname, Boolean is_na)

{
	if (blfsp1 == NULL)
		return FALSE;

	if (StringCmp(blfsp1->dbname, dbname) != 0)
		return FALSE; 
	if (blfsp1->is_prot == is_na)
		return FALSE; 
	if (blfsp1->bl3hp != bl3hp)
		return FALSE; 

	return TRUE;
}

/*********************************************************************

	Enables the fetching.  Initializes needed structures and calls
	BlastNetInit.

**********************************************************************/

Boolean LIBCALL 
BlastNetBioseqFetchEnable(BlastNet3Hptr bl3hp, CharPtr dbname, Boolean is_na, Boolean now)

{
        Boolean result;
        BlastNetFetchStructPtr blfsp = NULL;
        ObjMgrPtr omp;
        ObjMgrProcPtr ompp;

              /* check if already enabled ***/

        omp = ObjMgrGet();
        ompp = ObjMgrProcFind(omp, 0, "BlastNetBioseqFetch", OMPROC_FETCH);
        if (ompp != NULL)   /* already initialized */
	{
		blfsp = BlastNetFindFetchStruct(bl3hp, dbname, is_na);
			if (BlastNetFetchCompare(blfsp, bl3hp, dbname, is_na) == FALSE)
			{
				blfsp->dbname = MemFree(blfsp->dbname);
				blfsp->dbname = StringSave(dbname);
				blfsp->is_prot = (is_na == TRUE) ? FALSE : TRUE;
				blfsp->bl3hp = bl3hp;
				
			}
	}
	else
	{
		blfsp = BlastNetFindFetchStruct(bl3hp, dbname, is_na);

		ObjMgrProcLoad(OMPROC_FETCH, "BlastNetBioseqFetch", 
			"BlastNetBioseqFetch", OBJ_SEQID, 0,OBJ_BIOSEQ,0,
                        (Pointer)blfsp, BlastNetBioseqFetchFunc, PROC_PRIORITY_DEFAULT);

		blfsp->BlastNetFetchState = BLASTNET_INIT;
	}

	blfsp->ctr++;    /* count number of enables */

	if (blfsp->BlastNetFetchState == BLASTNET_READY)
	{
			  return TRUE;
	}

        if (now)
        {
		result = BlastNetInit(blfsp);
                if (! result)
                {
                        return result;
                }
		blfsp->BlastNetFetchState = BLASTNET_READY;
        }
        else
	{
		blfsp->BlastNetFetchState = BLASTNET_INIT;
	}

        return TRUE;
}

/*****************************************************************************
*
*		BlastNetBioseqFetchDisable()
*
*	Calls readdb_destruct if necessary to deallocate resources.
*
*****************************************************************************/
void LIBCALL BlastNetBioseqFetchDisable(BlastNet3Hptr bl3hp, CharPtr dbname, Boolean is_na)
{
        ObjMgrPtr omp;
        ObjMgrProcPtr ompp;
        BlastNetFetchStructPtr blfsp;

        omp = ObjMgrGet();
        ompp = ObjMgrProcFind(omp, 0, "BlastNetBioseqFetch", OMPROC_FETCH);
        if (ompp == NULL)   /* not initialized */
                return;

	blfsp = BlastNetFindFetchStruct(bl3hp, dbname, is_na);
	if (! blfsp->ctr)   /* no enables active */
		return;

	blfsp->ctr--;
	if (blfsp->ctr)   /* connection still pending */
			  return;

        if (blfsp->BlastNetFetchState == BLASTNET_READY)
	{
		blfsp->BlastNetFetchState = BLASTNET_DISABLE;  /* not active */
	}

        return;
}

/*
	Runs a BLAST request and returns a SeqAlignPtr for formatting.
	Note that the network connection must be established beforehand
	(i.e., BlastNet3BlockPtr blnet should be initialized).
*/
SeqAlignPtr LIBCALL
BlastBioseqNet(BlastNet3Hptr bl3hp, BioseqPtr bsp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback)

{
	return BlastBioseqNetCore(bl3hp, bsp, program, database, options, other_returns, error_returns, callback, NULL);
}

SeqAlignPtr LIBCALL
BlastSeqLocNet(BlastNet3Hptr bl3hp, SeqLocPtr slp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback)

{
	return BlastSeqLocNetCore(bl3hp, slp, program, database, options, other_returns, error_returns, callback, NULL);
}


SeqAlignPtr LIBCALL
BlastBioseqNetCore(BlastNet3Hptr bl3hp, BioseqPtr bsp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback, BLAST_MatrixPtr blast_matrix)

{
	BlastKABlkPtr ka_blk;
	BlastDbinfoPtr dbinfo;
	BlastMatrixPtr net_matrix;
	BLAST_MatrixPtr matrix;
	BlastNet3BlockPtr blnet;
	Boolean options_allocated = FALSE;
	CharPtr params_buffer;
	Int2 status;
	SeqAlignPtr seqalign = NULL;
	TxDfDbInfoPtr txdbinfo;
	ValNodePtr descr, mask;

	if (bl3hp == NULL || bsp == NULL || program == NULL || database == NULL)
		return NULL;

	if (error_returns)
		*error_returns = NULL;

	if (other_returns)
		*other_returns = NULL;

	/* If no options, use default. */
	if (options == NULL)
	{
		options = BLASTOptionNew(program, FALSE);
		options_allocated = TRUE;
	}

	status = BLASTOptionValidateEx(options, program, error_returns);
	if (status != 0)
	{	/* error messages in other_returns? */
		return NULL;
	}

	blnet = BlastNet3BlockNew(program, database);
	/* 
	Remove the Seq-descr as this is not needed for 
	BLASTing and the title often contains none-ASCII
	characters.  Keep the pointer and replace. 
	*/
	descr = bsp->descr;
	bsp->descr = NULL;

	blnet->bsp = bsp;
	blnet->parameters = BlastOptionsToParameters(options);
        if (options_allocated)
        {
                options = BLASTOptionDelete(options);
        }

	blnet->callback = callback;
	blnet->bl3hptr = bl3hp;
	blnet->blast_matrix = blast_matrix;

	seqalign = BlastBioseq(blnet, error_returns);
	if (other_returns)
	{
		*other_returns = NULL;
		mask = BlastGetMaskedLoc(blnet);
		if (mask)
			ValNodeLink(other_returns, mask);
		dbinfo = BlastGetDbInfo(blnet);
		txdbinfo = NetDbinfo2TxDbinfo(dbinfo);
		ValNodeAddPointer (other_returns, TXDBINFO, txdbinfo);
		dbinfo = BlastDbinfoFree(dbinfo);
		params_buffer = BlastGetParameterBuffer(blnet);
		ValNodeAddPointer(other_returns, TXPARAMETERS, params_buffer);
		ka_blk = BlastGetKaParams(blnet, FALSE);
		if (ka_blk)
			ValNodeAddPointer (other_returns, TXKABLK_NOGAP, ka_blk);
		ka_blk = BlastGetKaParams(blnet, TRUE);
		if (ka_blk)
			ValNodeAddPointer (other_returns, TXKABLK_GAP, ka_blk);
		net_matrix = NetBlastGetMatrix(blnet);
		matrix = BlastNetMatrixToBlastMatrix(net_matrix);
		if (matrix)
			ValNodeAddPointer (other_returns, TXMATRIX, matrix);
	}

	bsp->descr = descr;

	return seqalign;
}

SeqAlignPtr LIBCALL
BlastSeqLocNetCore(BlastNet3Hptr bl3hp, SeqLocPtr slp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback, BLAST_MatrixPtr blast_matrix)

{
	SeqAlignPtr	seqalign = NULL;
	ValNodePtr	vnp;

	SeqPortPtr	spp;
	Bioseq		bs;
	BioseqPtr	bsp = &bs;
	Int4		the_len, res_index;
	Int2		residue;
	ByteStorePtr	bp;
	SeqLocPtr	whole_slp;
	ValNodePtr	mask_loc;

	/* get Bioseq by SeqLoc */


	the_len = SeqLocLen(slp);
	bp = BSNew((Uint4)the_len);

	if (ISA_na(SeqLocMol(slp)))
	    spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
	else
	    spp = SeqPortNewByLoc (slp, Seq_code_iupacaa);

	MemSet((Pointer)bsp, 0, sizeof(Bioseq));
	bsp->length = (Int4)the_len;
	bsp->repr = Seq_repr_raw;
	bsp->mol = (Uint1) SeqLocMol(slp);
	bsp->id = SeqLocId(slp);

	if (ISA_na(bsp->mol))
	    bsp->seq_data_type = Seq_code_iupacna;
	else
	    bsp->seq_data_type = Seq_code_iupacaa;

	SeqPortSeek(spp, 0, SEEK_SET);
	BSSeek(bp, 0, SEEK_SET);
	for (res_index = 0; res_index < the_len; res_index++)
	{
	    residue = SeqPortGetResidue(spp);
	    BSPutByte(bp, residue);
	}
	bsp->seq_data = bp;

	SeqPortFree(spp);

	/* perform search */

	seqalign = BlastBioseqNet(bl3hp, bsp, program, database, options, 
		other_returns, error_returns, callback);
		
	/* offset the alignment */

	AdjustOffSetsInSeqAlign(seqalign, slp, NULL);
		
	if (other_returns) {
	    mask_loc = NULL;
	    for (vnp=*other_returns; vnp; vnp = vnp->next)
	    {
		switch (vnp->choice) {
		    case SEQLOC_MASKING_NOTSET:
		    case SEQLOC_MASKING_PLUS1:
		    case SEQLOC_MASKING_PLUS2:
		    case SEQLOC_MASKING_PLUS3:
		    case SEQLOC_MASKING_MINUS1:
		    case SEQLOC_MASKING_MINUS2:
		    case SEQLOC_MASKING_MINUS3:
			ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
			break;
		    default:
			break;
		}
	    }	

	    /* adjust offset in mask_loc */

	    if (mask_loc && slp) {
		SeqLocPtr	maskslp;
		SeqIntPtr	masksip;
		Int4	offset;
		ValNodePtr	vnp;

		whole_slp = NULL;

		ValNodeAddPointer(&whole_slp, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(bsp->id, SEQID_GI)));

		offset = GetOffsetInLoc(slp, whole_slp, SEQLOC_START);

		for (vnp = mask_loc; vnp; vnp = vnp->next) {

		    for (maskslp = (SeqLocPtr) vnp->data.ptrvalue; maskslp; maskslp = maskslp->next) {
			masksip = (SeqIntPtr) maskslp->data.ptrvalue;

			masksip->from += offset;
			masksip->to += offset;
		    }
		}
	    }
	}

	return seqalign;
}

#if EA
FILE *global_fp;
#endif

static  Boolean LIBCALLBACK
callback (BlastResponsePtr brp, Boolean PNTR cancel)

{

#if EA
        fprintf(global_fp, ".");
        fflush(global_fp);
#endif	
	return TRUE;
}

/*
	Formats a 'traditional' BLAST report.
*/

Boolean LIBCALL
TraditionalBlastReport(BioseqPtr bsp, BLAST_OptionsBlkPtr options, BlastNet3Hptr bl3hp, CharPtr program, CharPtr database, Boolean html, FILE *outfp, Boolean verbose, Uint4 print_options, Uint4 align_options, Int4 number_of_descriptions, Int4 number_of_alignments, Int4Ptr number_of_hits)

{
	BlastDbinfoPtr dbinfo;
	BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
	BlastPruneSapStructPtr prune;
	BLAST_MatrixPtr matrix;
	Boolean query_is_na, db_is_na;
	CharPtr params_buffer=NULL;
	Int4 number_of_hits_private=0;
	SeqAlignPtr seqalign;
        SeqAnnotPtr seqannot=NULL;
	TxDfDbInfoPtr tx_dbinfo=NULL, tx_dbinfo_head;
	ValNodePtr mask_loc, other_returns, error_returns, vnp;
        Uint1 align_type;
        Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];

        MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
        MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));

	if (bsp == NULL || bl3hp == NULL || program == NULL || database == NULL || outfp == NULL)
		return FALSE;

#if EA
NlmMutexLockEx(&formating_mutex);
#endif

	align_type = BlastGetTypes(program, &query_is_na, &db_is_na);
#if EA
	global_fp = outfp;
#endif

       	init_buff_ex(85);
	dbinfo = BlastRequestDbInfo(bl3hp, database, !db_is_na);
	if (dbinfo)
       		PrintDbInformationBasic(database, !db_is_na, 70, dbinfo->definition, dbinfo->number_seqs, dbinfo->total_length, outfp, html);
	dbinfo = BlastDbinfoFree(dbinfo);
       	free_buff();
#if EA
NlmMutexUnlock(formating_mutex);
#endif

	fprintf(outfp, "Searching");

	seqalign = BlastBioseqNet(bl3hp, bsp, program, database, options, &other_returns, &error_returns, callback);

	fprintf(outfp, "done");
		
	BlastErrorPrintExtra(error_returns, TRUE, outfp);
	
	mask_loc = NULL;
	for (vnp=other_returns; vnp; vnp = vnp->next)
	{
			switch (vnp->choice) {
				case TXDBINFO:
					tx_dbinfo = vnp->data.ptrvalue;
					break;
				case TXKABLK_NOGAP:
					ka_params = vnp->data.ptrvalue;
					break;
				case TXKABLK_GAP:
					ka_params_gap = vnp->data.ptrvalue;
					break;
				case TXPARAMETERS:
					params_buffer = vnp->data.ptrvalue;
					break;
				case TXMATRIX:
					matrix = vnp->data.ptrvalue;
					break;
				case SEQLOC_MASKING_NOTSET:
				case SEQLOC_MASKING_PLUS1:
				case SEQLOC_MASKING_PLUS2:
				case SEQLOC_MASKING_PLUS3:
				case SEQLOC_MASKING_MINUS1:
				case SEQLOC_MASKING_MINUS2:
				case SEQLOC_MASKING_MINUS3:
					ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
					break;
				default:
					break;
			}
	}	

#if EA
NlmMutexLockEx(&formating_mutex);
#endif

	if (seqalign)
	{
		seqannot = SeqAnnotNew();
        	seqannot->type = 2;
		AddAlignInfoToSeqAnnot(seqannot, align_type);
        	seqannot->data = seqalign;
		prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_descriptions, NULL);
		ObjMgrSetHold();
       		init_buff_ex(85);
                PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, print_options, FIRST_PASS, NULL);
                free_buff();

		prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_alignments, prune);
		seqannot->data = prune->sap;
		if (align_options & TXALIGN_MASTER)
			ShowTextAlignFromAnnot(seqannot, 60, outfp, f_order, g_order, align_options, NULL, mask_loc, NULL);
		else
			ShowTextAlignFromAnnot(seqannot, 60, outfp, f_order, g_order, align_options, NULL, mask_loc, FormatScoreFunc);
		seqannot->data = seqalign;
		number_of_hits_private = prune->original_number; 
		prune = BlastPruneSapStructDestruct(prune);
		ObjMgrClearHold();
	}

	if (verbose)
	{
       		init_buff_ex(85);
		tx_dbinfo_head = tx_dbinfo;
		while (tx_dbinfo)
		{
                	PrintDbReport(tx_dbinfo, 70, outfp);
			tx_dbinfo = tx_dbinfo->next;
		}
		tx_dbinfo_head = TxDfDbInfoDestruct(tx_dbinfo_head);

		if (ka_params)
		{
                	PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
			MemFree(ka_params);
		}

		if (ka_params_gap)
		{
                	PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
			MemFree(ka_params_gap);
		}

                PrintTildeSepLines(params_buffer, 70, outfp);
                MemFree(params_buffer);
                free_buff();
	}

#if EA
NlmMutexUnlock(formating_mutex);
#endif
	if (seqannot)
		seqannot = SeqAnnotFree(seqannot);

	if (number_of_hits)
		*number_of_hits = number_of_hits_private;

	return TRUE;
}

Boolean LIBCALL
TraditionalBlastReportLoc(SeqLocPtr slp, BLAST_OptionsBlkPtr options, BlastNet3Hptr bl3hp, CharPtr program, CharPtr database, Boolean html, FILE *outfp, Boolean verbose, Uint4 print_options, Uint4 align_options, Int4 number_of_descriptions, Int4 number_of_alignments, Int4Ptr number_of_hits)

{
	BlastDbinfoPtr dbinfo;
	BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
	BlastPruneSapStructPtr prune;
	BLAST_MatrixPtr matrix;
	Boolean query_is_na, db_is_na;
	CharPtr params_buffer=NULL;
	SeqAlignPtr seqalign;
        SeqAnnotPtr seqannot=NULL;
	TxDfDbInfoPtr tx_dbinfo=NULL, tx_dbinfo_head;
	ValNodePtr mask_loc, other_returns, error_returns, vnp;
        Uint1 align_type;
        Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];

        MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
        MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));


	if (slp == NULL || bl3hp == NULL || program == NULL || database == NULL || outfp == NULL)
		return FALSE;

	align_type = BlastGetTypes(program, &query_is_na, &db_is_na);
#if EA
	global_fp = outfp;
#endif

       	init_buff_ex(85);
	dbinfo = BlastRequestDbInfo(bl3hp, database, !db_is_na);
	if (dbinfo)
       		PrintDbInformationBasic(database, !db_is_na, 70, dbinfo->definition, dbinfo->number_seqs, dbinfo->total_length, outfp, html);
	dbinfo = BlastDbinfoFree(dbinfo);
       	free_buff();

	fprintf(outfp, "Searching");

	seqalign = BlastSeqLocNet(bl3hp, slp, program, database, options, &other_returns, &error_returns, callback);

	fprintf(outfp, "done");

	mask_loc = NULL;

	BlastErrorPrintExtra(error_returns, TRUE, outfp);

	for (vnp=other_returns; vnp; vnp = vnp->next)
	{
			switch (vnp->choice) {
				case TXDBINFO:
					tx_dbinfo = vnp->data.ptrvalue;
					break;
				case TXKABLK_NOGAP:
					ka_params = vnp->data.ptrvalue;
					break;
				case TXKABLK_GAP:
					ka_params_gap = vnp->data.ptrvalue;
					break;
				case TXPARAMETERS:
					params_buffer = vnp->data.ptrvalue;
					break;
				case TXMATRIX:
					matrix = vnp->data.ptrvalue;
					break;
				case SEQLOC_MASKING_NOTSET:
				case SEQLOC_MASKING_PLUS1:
				case SEQLOC_MASKING_PLUS2:
				case SEQLOC_MASKING_PLUS3:
				case SEQLOC_MASKING_MINUS1:
				case SEQLOC_MASKING_MINUS2:
				case SEQLOC_MASKING_MINUS3:
#if 1					
					ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
#endif					
					break;
				default:
					break;
			}
	}	

	if (seqalign)
	{

		seqannot = SeqAnnotNew();
        	seqannot->type = 2;
		AddAlignInfoToSeqAnnot(seqannot, align_type);
        	seqannot->data = seqalign;
		prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_descriptions, NULL);
		ObjMgrSetHold();
       		init_buff_ex(85);
                PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, print_options, FIRST_PASS, NULL);
                free_buff();

		prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_alignments, prune);
		seqannot->data = prune->sap;

		if (align_options & TXALIGN_MASTER)
			ShowTextAlignFromAnnot(seqannot, 60, outfp, f_order, g_order, align_options, NULL, mask_loc, NULL);
		else
			ShowTextAlignFromAnnot(seqannot, 60, outfp, f_order, g_order, align_options, NULL, mask_loc, FormatScoreFunc);
		seqannot->data = seqalign;
		prune = BlastPruneSapStructDestruct(prune);
		ObjMgrClearHold();
	}

	if (verbose)
	{
       		init_buff_ex(85);
		tx_dbinfo_head = tx_dbinfo;
		while (tx_dbinfo)
		{
                	PrintDbReport(tx_dbinfo, 70, outfp);
			tx_dbinfo = tx_dbinfo->next;
		}
		tx_dbinfo_head = TxDfDbInfoDestruct(tx_dbinfo_head);

		if (ka_params)
		{
                	PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
			MemFree(ka_params);
		}

		if (ka_params_gap)
		{
                	PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
			MemFree(ka_params_gap);
		}

                PrintTildeSepLines(params_buffer, 70, outfp);
                MemFree(params_buffer);
                free_buff();
	}

	if (seqannot)
		seqannot = SeqAnnotFree(seqannot);


	return TRUE;
}


/*
        Converst the BlastParametersPtr (used by network service) to
        BLAST_OptionsBlkPtr (used by blast).
*/

BLAST_OptionsBlkPtr
parametersToOptions (BlastParametersPtr parameters, CharPtr program, ValNodePtr PNTR error_return)

{
        BLAST_OptionsBlkPtr options;
        Int2 status;

        if (program == NULL)
                return NULL;

        if (parameters == NULL)
        {
                options = BLASTOptionNew(program, TRUE);
        }
        else
        {
                options = BLASTOptionNew(program, (Boolean) parameters->gapped_alignment);
                options->threshold_first = parameters->first_threshold;
                options->threshold_second = parameters->second_threshold;
                if (parameters->Cutoff_cutoff)
                {
                        if (parameters->Cutoff_cutoff->choice == Cutoff_cutoff_evalue)
                                options->expect_value = parameters->Cutoff_cutoff->data.realvalue;
                        else if (parameters->Cutoff_cutoff->choice == Cutoff_cutoff_score)
                                options->cutoff_s = parameters->Cutoff_cutoff->data.intvalue;
                }
                if (parameters->Cutoff2_cutoff2)
                {
                        if (parameters->Cutoff2_cutoff2->choice == Cutoff2_cutoff2_evalue)
                                options->e2 = parameters->Cutoff2_cutoff2->data.realvalue;
                        else if (parameters->Cutoff2_cutoff2->choice == Cutoff2_cutoff2_score)
                                options->cutoff_s2 = parameters->Cutoff2_cutoff2->data.intvalue;
                }
                /* compensates for client not providing this.  Remove this at some point? */
                if (parameters->hitlist_size != 0)
                        options->hitlist_size = parameters->hitlist_size;
                options->penalty = parameters->nucl_penalty;
                options->reward = parameters->nucl_reward;
                options->gap_open = parameters->gap_open;
                options->gap_extend = parameters->gap_extend;
                /* compensates for client not providing this.  Remove this at some point? */
                if (parameters->genetic_code != 0)
                        options->genetic_code = parameters->genetic_code;
                /* compensates for client not providing this.  Remove this at some point? */
                if (parameters->db_genetic_code != 0)
                        options->db_genetic_code = parameters->db_genetic_code;

                options->filter = parameters->low_complexity_filtering;
                options->ethresh = parameters->ethresh;
                options->maxNumPasses = parameters->max_num_passes;
                options->pseudoCountConst = parameters->pseudo_count_const;
                
                options->gifile = StringSave(parameters->gifile);
                options->gilist = parameters->gilist;
		if (parameters->matrix)
                	options->matrix = StringSave(parameters->matrix);
		if (parameters->filter_string)
               		options->filter_string = StringSave(parameters->filter_string);
                if (parameters->entrez_query)
                    options->entrez_query = StringSave(parameters->entrez_query);
                /* compensates for client not providing this.  Remove this at some point? */
		if (parameters->word_size)
                    options->wordsize = parameters->word_size;
        }

	if (status = BLASTOptionValidateEx(options, program, error_return)) {
            return NULL;
	}

        return options;
}