/*   tax0.c
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
* File Name:  tax0.c
*
* Author:  Soussov
*
* Version Creation Date:   01/22/96
*
* $Revision: 6.0 $
*
* File Description: 
*       API for Taxonomy Archive service
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*
*
* RCS Modification History:
* $Log: tax0.c,v $
* Revision 6.0  1997/08/25 18:41:19  madden
* Revision changed to 6.0
*
* Revision 5.2  1996/07/03 19:02:47  kans
* removed inclusion of taxerr.h
*
 * Revision 5.1  1996/07/03  18:46:16  epstein
 * eliminate double-quotes in favor of angle-brackets and eliminate mappings.h
 *
 * Revision 5.0  1996/05/28  14:15:13  ostell
 * Set to revision 5.0
 *
 * Revision 1.1  1996/03/06  17:05:19  soussov
 * Initial revision
 *
*/

/** for ErrPostEx() ****/

static char *this_module = "tax0";

#ifdef THIS_MODULE
#undef THIS_MODULE
#endif

#define THIS_MODULE this_module
static char *this_file = __FILE__;
#define THIS_FILE this_file

#include <ncbinet.h>
#include <objall.h>
#include <objfeat.h>
#include <tax0.h>
#include <taxinc.h>

static Taxon0RespPtr NetTaxArchReadAsn PROTO((void));
static Boolean ReestablishNetTaxArch PROTO((void));
static Boolean NetInit PROTO((void));
static Boolean ForceNetInit PROTO((void));
static Boolean NetFini PROTO((void));
static Boolean GenericReestablishNet PROTO((CharPtr svcName, Boolean showErrs));

static TaxonIdListPtr s_TaxArchGetTaxId PROTO((TaxonNamePtr tnp));
static Pointer s_TaxArchAccessServer PROTO((Int4 id_tax, Int2 request, Int2 resp));
static Pointer s_TaxArchGetFromServer PROTO((Int4 id_tax, Int2 req, Int2 resp));

static NI_HandPtr svcp = NULL;
static AsnIoPtr   asnin = NULL;
static AsnIoPtr   asnout = NULL;
Taxon0ReqPtr taxrp;
Taxon0RespPtr taxbp;
static Boolean num_attached = 0;
static Boolean reallyFinal = TRUE;
static NI_DispatcherPtr dispatcher;
static Boolean (*myNetInit) PROTO((void));

/*****************************************************************************
*
*   Tax0Init ()
*
*****************************************************************************/

Boolean Tax0Init (void)

{
    DataVal    av;

    myNetInit = Tax0Init;

    if (! NetInit())
        return FALSE;

    svcp = NI_GenericGetService(dispatcher, NULL, "TAXONOMY", "Taxonomy", TRUE);
    if (svcp == NULL)
    {
        ErrPostEx(SEV_ERROR, 0, 0, "NI_ServiceGet [%s] (%s)", ni_errlist[ni_errno], ni_errtext);
        Tax0Fini();
        return FALSE;
    }

    asnin = svcp->raip;
    asnout = svcp->waip;

    /**********************************************************/

    taxrp = ValNodeNew(NULL);
    taxrp->choice = Taxon0Req_init;
    Taxon0ReqAsnWrite (taxrp, asnout, NULL);
    AsnIoReset(asnout);
    Taxon0ReqFree (taxrp);

    if ((taxbp = NetTaxArchReadAsn()) == NULL)
    {
        return FALSE;
    }
    else
    {
        taxbp->data.ptrvalue = NULL;
        Taxon0RespFree (taxbp);
        return TRUE;
    }
}

/*****************************************************************************
*
*   Tax0Fini ()
*
*****************************************************************************/

static Boolean s_TaxArchFini (void)

{
    Boolean retval = TRUE;

    if (asnout != NULL && asnin != NULL)
    {
        taxrp = ValNodeNew(NULL);
        taxrp->choice = Taxon0Req_fini;
        Taxon0ReqAsnWrite (taxrp, asnout, NULL);
        AsnIoReset(asnout);
        Taxon0ReqFree (taxrp);

        if ((taxbp = NetTaxArchReadAsn()) == NULL)
        {
            retval = FALSE;
        }
        else
        {
            taxbp->data.ptrvalue = NULL;
            Taxon0RespFree (taxbp);
        }
    }

    NetFini();

    return retval;
}

/* the only thing done here is to suppress errors */

Boolean Tax0Fini (void)

{
    short erract;
    ErrDesc err;
    Boolean retval;

    ErrGetOpts(&erract, NULL);
    ErrSetOpts(ERR_IGNORE, 0);
    ErrFetch(&err);

    retval = s_TaxArchFini();

    ErrSetOpts(erract, 0);
    ErrFetch(&err);

    return retval;
}


/*****************************************************************************
*
*   s_TaxArchGetTaxComplete ()
*
*****************************************************************************/
static
TaxCompleteListPtr s_TaxArchGetTaxComplete( TaxonIdNamePtr tinp )
{
    Taxon0ReqPtr taxrp;
    TaxCompleteListPtr tclp = NULL;

    taxrp = ValNodeNew( NULL );
    taxrp->choice = Taxon0Req_getcomplete;
    taxrp->data.ptrvalue = (Pointer) tinp;
    Taxon0ReqAsnWrite( taxrp, asnout, NULL );
    AsnIoReset( asnout );
    taxrp->data.ptrvalue = NULL;
    Taxon0ReqFree ( taxrp );

    if ( ( taxbp = NetTaxArchReadAsn() ) == NULL )
        return NULL;

    if ( taxbp->choice != Taxon0Resp_getcomplete )
    {
        Taxon0RespFree( taxbp );
        return NULL;
    }
    tclp = (TaxCompleteListPtr) (taxbp->data.ptrvalue);
    taxbp->data.ptrvalue = NULL;
    Taxon0RespFree( taxbp );

    return tclp;
}


/*****************************************************************************
*
*   s_TaxArchGetTaxId ()
*
*****************************************************************************/
static
TaxonIdListPtr s_TaxArchGetTaxId( TaxonNamePtr tnp )
{
    Taxon0ReqPtr taxrp;
    TaxonIdListPtr tilp = NULL;

    taxrp = ValNodeNew( NULL );
    taxrp->choice = Taxon0Req_getid;
    taxrp->data.ptrvalue = (Pointer) tnp;
    Taxon0ReqAsnWrite( taxrp, asnout, NULL );
    AsnIoReset( asnout );
    taxrp->data.ptrvalue = NULL;
    Taxon0ReqFree ( taxrp );

    if ( ( taxbp = NetTaxArchReadAsn() ) == NULL )
        return NULL;

    if ( taxbp->choice != Taxon0Resp_getid )
    {
        Taxon0RespFree( taxbp );
        return NULL;
    }
    tilp = (TaxonIdListPtr) (taxbp->data.ptrvalue);
    taxbp->data.ptrvalue = NULL;
    Taxon0RespFree( taxbp );

    return tilp;
}


TaxonIdListPtr Tax0GetTaxId( TaxonNamePtr tnp )
{
    Int4 i;
    short erract;
    ErrDesc err;
    TaxonIdListPtr tilp = NULL;

    for (i = 0; i < TAXARCH_SERV_RETRIES; i++)
    {
        if (i > 0)
        {
            if (! ReestablishNetTaxArch())
                break;
        }
    
        ErrGetOpts(&erract, NULL);
        ErrSetOpts(ERR_IGNORE, 0);
        ErrFetch(&err);

        tilp = s_TaxArchGetTaxId( tnp );

        ErrSetOpts(erract, 0);
        if (! ErrFetch(&err))
            break; /* success */
    }

    return( tilp );
}

/*****************************************************************************
*
*  s_TaxArchAccessServer ()
*
*****************************************************************************/
static
Pointer s_TaxArchAccessServer( Int4 id_tax, Int2 request, Int2 resp )
{
    Taxon0ReqPtr taxrp;
    Pointer retptr = NULL;

    taxrp = ValNodeNew( NULL );
    taxrp->choice = request;
    taxrp->data.intvalue = id_tax;
    Taxon0ReqAsnWrite( taxrp, asnout, NULL );
    AsnIoReset( asnout );
    Taxon0ReqFree ( taxrp );

    if ( ( taxbp = NetTaxArchReadAsn() ) == NULL )
        return NULL;

    if ( taxbp->choice != resp )
    {
        Taxon0RespFree( taxbp );
        return NULL;
    }

    retptr = taxbp->data.ptrvalue;
    taxbp->data.ptrvalue = NULL; /* for clean free */
    Taxon0RespFree( taxbp );

    return retptr;
}

/*****************************************************************************
*
*   s_TaxArchGetFromServer()
*
*****************************************************************************/
static
Pointer s_TaxArchGetFromServer( Int4 id_tax, Int2 req, Int2 resp )
{
    Int4 i;
    short erract;
    ErrDesc err;
    Pointer retptr = NULL;

    for (i = 0; i < TAXARCH_SERV_RETRIES; i++)
    {
        if (i > 0)
        {
            if (! ReestablishNetTaxArch())
                break;
        }
    
        ErrGetOpts( &erract, NULL );
        ErrSetOpts( ERR_IGNORE, 0 );
        ErrFetch( &err );

        retptr = s_TaxArchAccessServer( id_tax, req, resp );

        ErrSetOpts( erract, 0 );
        if (! ErrFetch( &err ) )
            break; /* success */
    }

    return( retptr );
}

/*****************************************************************************
*
*   Tax0GetChildren ()
*
*****************************************************************************/

TaxonIdListPtr Tax0GetChildren( Int4 id_tax )
{
    TaxonIdListPtr tilp = NULL;

    tilp = (TaxonIdListPtr) s_TaxArchGetFromServer( id_tax,
                                  Taxon0Req_getchildren, Taxon0Resp_gettaxon );
    return( tilp );
}

/*****************************************************************************
*
*   Tax0GetParents ()
*
*****************************************************************************/

TaxonIdListPtr Tax0GetParents( Int4 id_tax )
{
    TaxonIdListPtr tilp = NULL;

    tilp = (TaxonIdListPtr) s_TaxArchGetFromServer( id_tax,
                                  Taxon0Req_getparents, Taxon0Resp_gettaxon );
    return( tilp );
}


/*****************************************************************************
*
*   Tax0GetRef ()
*
*****************************************************************************/

OrgRefPtr Tax0GetRef( Int4 id_tax )
{
    OrgRefPtr orp = NULL;

    orp = (OrgRefPtr) s_TaxArchGetFromServer( id_tax, 
                              Taxon0Req_getref, Taxon0Resp_getref );
    return( orp );
}


/*****************************************************************************
*
*   Tax0GetComplete ()
*
*****************************************************************************/

TaxCompleteListPtr Tax0GetComplete( TaxonIdNamePtr tinp )
{
    TaxCompleteListPtr tclp = NULL;
    TaxCompletePtr tcp;
    int i;

    tclp = (TaxCompleteListPtr) s_TaxArchGetTaxComplete( tinp );

    if (tclp == NULL) {
      return NULL;
    }

    for ( i = 0, tcp = tclp->info; i < tclp->num; tcp = tcp->next, i++ ) {
        if ( tcp->comname[0]  == '\0') tcp->comname   = MemFree(tcp->comname);
        if ( tcp->synonyms[0] == '\0') tcp->synonyms  = MemFree(tcp->synonyms);
        if ( tcp->name_gc[0]  == '\0') tcp->name_gc   = MemFree(tcp->name_gc);
        if ( tcp->name_mgc[0] == '\0') tcp->name_mgc  = MemFree(tcp->name_mgc);
        if ( tcp->gb_div[0]   == '\0') tcp->gb_div    = MemFree(tcp->gb_div);
        if ( tcp->embl_code[0]== '\0') tcp->embl_code = MemFree(tcp->embl_code);
    }
    return( tclp );
}

/*****************************************************************************
*
*   NetTaxArchReadAsn ()
*
*****************************************************************************/

static Taxon0RespPtr NetTaxArchReadAsn(void)
{
    Taxon0RespPtr taxbp;
    short erract;
    ErrDesc err;

    ErrGetOpts(&erract, NULL);
    ErrSetOpts(ERR_CONTINUE, 0);
    ErrFetch(&err); /* clear any pending error */

    taxbp = Taxon0RespAsnRead(asnin, NULL);

    if (ErrFetch(&err))
    {
        ErrPost (CTX_UNKNOWN, 1, "Null message read from server");
    }
    ErrSetOpts(erract, 0);

    return taxbp;
}

/*****************************************************************************
*
*   ReestablishNetTaxArch ()
*
*****************************************************************************/

static Boolean ReestablishNetTaxArch(void)
{
    return GenericReestablishNet("Taxonomy", TRUE);
}

/*****************************************************************************
*
*   GenericReestablishNet ()
*
*****************************************************************************/

static Boolean GenericReestablishNet(CharPtr svcName, Boolean showErrs)
{
    MonitorPtr mon = NULL;
    Boolean retval;
    CharPtr buf;

    buf = MemNew(2 * StrLen(svcName) + 60);

    if (showErrs) {
        sprintf (buf, "Re-establishing %s Service", svcName);
        mon = MonitorStrNew(buf, 40);
        sprintf (buf, "Requesting %s service", svcName);
        MonitorStrValue(mon, buf);
    }
    NetFini();
    retval = TRUE;

    if (! myNetInit())
    {
        sprintf (buf, "%s get failed; re-contacting dispatcher", svcName);
        MonitorStrValue(mon, buf);
        retval = FALSE;
        if (ForceNetInit())
        { /* successfully established contact w/dispatcher */
            sprintf (buf, "%s get failed; re-requesting %s service",
                     svcName, svcName);
            MonitorStrValue(mon, buf);
            retval = myNetInit();
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
*   NetInit ()
*
*****************************************************************************/

static Boolean
NetInit(void)
{
    if (num_attached++ > 0)
        return TRUE;

    return ((dispatcher = NI_GenericInit(NULL, NULL, TRUE, NULL, 0)) != NULL);
}


/*****************************************************************************
*
*   ForceNetInit ()
*
*****************************************************************************/

static Boolean ForceNetInit(void)
{
    Boolean retval;

    reallyFinal = FALSE;
    num_attached = 0; /* force re-attempt to contact dispatcher */
    retval = NetInit();
    reallyFinal = TRUE;

    return retval;
}

/*****************************************************************************
*
*   NetFini ()
*
*****************************************************************************/

static Boolean NetFini(void)
{
   if (num_attached > 0)
        num_attached--;

    if (num_attached == 0)
    {
        NI_ServiceDisconnect(svcp);
        svcp = NULL;
        NI_EndServices (dispatcher);
        dispatcher = NULL;
    }

    return TRUE;
}


OrgRefPtr TaxArchGetRef(Int4 id_tax)
{
  Taxon1DataPtr res;
  OrgRefPtr orp;

  res= tax1_getbyid(id_tax);
  if(res == NULL) return NULL;
  orp= res->org;
  res->org= NULL;
  Taxon1DataFree(res);
  return orp;
}

