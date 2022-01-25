/*
*
*
* RCS Modification History:
* $Log: taxsmain.c,v $
* Revision 6.1  1998/02/19 20:03:04  soussov
* -z argument added
*
* Revision 6.0  1997/08/25 18:42:00  madden
* Revision changed to 6.0
*
* Revision 1.1  1996/05/23 16:55:55  soussov
* Initial revision
*
 * Revision 1.7  1995/05/17  17:58:52  epstein
 * add RCS log revision history
 *
*/

#include <varargs.h>
#include <ncbi.h>
#include <ncbinet.h>
#include <accentr.h>
#include <objtaxon.h>
#include <tax_cmmn.h>

TaxonIdListPtr FsTaxGetTaxId( TaxonNamePtr tnp );
TaxonIdListPtr FsTaxGetParents( int id_tax );
TaxonIdListPtr FsTaxGetChildren( int id_tax );
OrgRefPtr FsTaxGetRef( int id_tax );
ValNodePtr FsTaxGetTaxonline( int id_tax );
ValNodePtr FsTaxGetWithinDiv( int id_tax );
GeneticCodeListPtr FsTaxGetGeneticCode( int id_tax );
TaxCompleteListPtr FsTaxGetComplete( TaxonIdNamePtr tinp );


TaxSendResp( Int2 resp_choice, Pointer ptr, AsnIoPtr asnout )
{
    Taxon0RespPtr taxbp;

    taxbp = ValNodeNew(NULL);
    if ( ptr == NULL )
    {
	taxbp->choice = Taxon0Resp_error;
	taxbp->data.intvalue = 0;
    }
    else
    {
    	taxbp->choice = resp_choice;
    	taxbp->data.ptrvalue = ptr;
    }
    Taxon0RespAsnWrite( taxbp, asnout, NULL );

    taxbp->data.ptrvalue = NULL;
    Taxon0RespFree( taxbp );
}

#define VBUF_SIZ 1024

#ifdef DO_LOGGING
void 
LOGIT(int va_alist)
{
    va_list             args;
    CharPtr             message;
    Char                msgstr[VBUF_SIZ];
    time_t              now;
    Char                ctimestr[30];
    
    MemFill(msgstr, NULL, VBUF_SIZ);
    va_start(args);
    message = va_arg(args, CharPtr);
    if (message != NULL)
        vsprintf(msgstr, message, args);
    va_end(args);
    time (&now);
    StrCpy (ctimestr, ctime(&now));
    ctimestr[StrLen(ctimestr) - 1] = '\0'; /* remove trailing newline */
    fprintf(stderr, "%s : %s\n", ctimestr,msgstr);
}
#else
void 
LOGIT(int va_alist)
{
    va_list             args;
    CharPtr             message;
    Char                msgstr[VBUF_SIZ];
    time_t              now;
    Char                ctimestr[30];
    
    va_start(args);
    va_end(args);
}
#endif



main(int argc, char *argv[])
{
    Taxon0ReqPtr taxrp;
    Taxon0RespPtr taxbp;
    Boolean debug = FALSE;
    Boolean filter= FALSE;
    short erract;
    ErrDesc err;
    int arg;
    NI_HandPtr hp;
    Char buf[100];
    AsnIoPtr asnin;
    AsnIoPtr asnout;
    int read_timeout;
    Boolean done = FALSE;

    TaxonIdListPtr tilp;
    OrgRefPtr orp;
    ValNodePtr lineage;
    ValNodePtr within_div;
    GeneticCodeListPtr gclp;
    TaxCompleteListPtr tclp;
    TaxonIdNamePtr tinp;
    TaxonNamePtr tnp;

		ErrSetFatalLevel(SEV_MAX);

    if (argc > 1)
    {
        arg = 1;
        if (StrCmp(argv[1], "-d") == 0)
        {
            arg++;
            debug = TRUE;
        }
	else if(StrCmp(argv[1], "-z") == 0) {
	    arg++;
	    filter= TRUE;
	}
    }

    if (getenv("NI_LOCAL_ACCESS_DENIED") != NULL)
    {
        if (!debug)
            NI_ServerNACK("Unauthorized access attempt for Taxonomy server");
        return -1;
    }

    LOGIT("Entering InitTaxDB");
    if (! InitTaxDB() )
    {   
        LOGIT("InitTaxDB failed");

        if (!debug)
            NI_ServerNACK("Taxonomy service: Unable to initialize DB");
        return( -1 );
    }
    LOGIT("InitTaxDB OK");


    if (!debug) {
	if(filter) {
	    asnin = AsnIoNew(ASNIO_BIN_IN, stdin, NULL, NULL, NULL);
	    asnout = AsnIoNew(ASNIO_BIN_OUT, stdout, NULL, NULL, NULL);
	}
	else {
	    NI_ServerACK();
    
	    hp = NI_OpenASNIO();

	    /* this read-timeout is effectively an idle timeout for */
	    /* the server process; the process will terminate upon  */
	    /* read-timeout                                         */
	    GetAppParam("NCBI", "NET_SERV", "SERV_INACT_TIMER", "10",
			buf, sizeof buf);
	    read_timeout = atoi(buf) * 60; /* param is minutes */
	    MsgSetReadTimeout(hp, read_timeout);

	    asnin = hp->raip;
	    asnout = hp->waip;
	}
    } else {
	asnin = AsnIoOpen("taxserv.inp", "r");
        asnout = AsnIoOpen("taxserv.out", "w");
    }

    while (!done) {
        /* encountering EOF on reading is a "normal" occurrence, */
        /* and does not merit an error message                   */
        ErrGetOpts(&erract, NULL);
        ErrSetOpts(ERR_IGNORE, 0);
        ErrFetch(&err); /* clear any pending error, which can be ignored */

        taxrp = Taxon0ReqAsnRead(asnin, NULL);

        if (ErrFetch(&err))
        {
	    done = TRUE;
            ErrPostEx(SEV_ERROR,1,1, "Error encountered on AsnReadId %d", err);
            break; /* client terminated */
        }
        ErrSetOpts(erract, 0);

        if (taxrp == NULL)
        {
	    done = TRUE;
            ErrPostEx(SEV_ERROR,1,1, "Null AsnReadId");
            break; /* client terminated */
        }

        switch (taxrp->choice) {

        case Taxon0Req_init:

            taxbp = ValNodeNew(NULL);
            taxbp->choice = Taxon0Resp_init;
            taxbp->data.ptrvalue = NULL;
            Taxon0RespAsnWrite (taxbp, asnout, NULL);
            Taxon0RespFree (taxbp);
            taxrp->data.ptrvalue = NULL;
            break;

        case Taxon0Req_getid:

	    tnp = (TaxonNamePtr) taxrp->data.ptrvalue;
	    tilp = FsTaxGetTaxId( tnp );

	    taxrp->data.ptrvalue = NULL;

	    TaxonNameFree( tnp );
	    TaxSendResp( Taxon0Resp_getid, (Pointer) tilp, asnout );

            break;

        case Taxon0Req_getref:

	    orp = FsTaxGetRef( taxrp->data.intvalue );

	    TaxSendResp( Taxon0Resp_getref, (Pointer) orp, asnout );

	    if ( orp != NULL && orp->syn != NULL )
		orp->syn = ValNodeFree( orp->syn );

            break;

        case Taxon0Req_getchildren:

	    tilp = FsTaxGetChildren( taxrp->data.intvalue );

	    TaxSendResp( Taxon0Resp_gettaxon, (Pointer) tilp, asnout );

            break;

        case Taxon0Req_getparents:

            tilp = FsTaxGetParents( taxrp->data.intvalue );

            TaxSendResp( Taxon0Resp_gettaxon, (Pointer) tilp, asnout );

            break;

        case Taxon0Req_getgeneticcode:

	    gclp = FsTaxGetGeneticCode( taxrp->data.intvalue );

            TaxSendResp( Taxon0Resp_getgeneticcode, (Pointer) gclp, asnout );

            break;

        case Taxon0Req_gettaxonline:

	    lineage = FsTaxGetTaxonline( taxrp->data.intvalue );

            TaxSendResp( Taxon0Resp_gettaxonline, 
                            (Pointer)(lineage?lineage->data.ptrvalue:""), asnout );

            break;

        case Taxon0Req_getdivision:

	    within_div = FsTaxGetWithinDiv( taxrp->data.intvalue );

            TaxSendResp( Taxon0Resp_getdivision, 
                            (Pointer) (within_div?within_div->data.ptrvalue:""), asnout );

            break;

        case Taxon0Req_getcomplete:

	    tinp = (TaxonIdNamePtr)taxrp->data.ptrvalue;
	    tclp = FsTaxGetComplete( tinp );

	    taxrp->data.ptrvalue = NULL;

            TaxSendResp( Taxon0Resp_getcomplete, (Pointer) tclp, asnout );

	    TaxonIdNameFree( tinp );
	    TaxCompleteListFree( tclp );

            break;

        case Taxon0Req_fini:

            done = TRUE;
            taxbp = ValNodeNew(NULL);
            taxbp->choice = Taxon0Resp_fini;
            taxbp->data.ptrvalue = NULL;

            Taxon0RespAsnWrite (taxbp, asnout, NULL);
            Taxon0RespFree (taxbp);
            taxrp->data.ptrvalue = NULL;
            break;
        }

        AsnIoReset (asnout);
        Taxon0ReqFree (taxrp);
    }

    CloseTaxDB();

    AsnIoClose (asnin);
    AsnIoClose (asnout);
}
