/*   tax3api.c
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
* File Name:  tax3api.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   7/8/04
*
* $Revision: 1.35 $
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

#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <tax3api.h>
#include <sqnutils.h>
#include <subutil.h>
#include <findrepl.h>

/* low-level connection functions */

static Boolean text_tax_asn = FALSE;
static Boolean text_tax_set = FALSE;

NLM_EXTERN CONN Tax3OpenConnection (
  void
)

{
#ifdef OS_UNIX
  CharPtr  str;

  if (! text_tax_set) {
    str = (CharPtr) getenv ("TEXT_TAX_ASN");
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "TRUE") == 0) {
        text_tax_asn = TRUE;
      }
    }
    text_tax_set = TRUE;
  }
#endif

  return QUERY_OpenServiceQuery (text_tax_asn ? "TaxService3Text" : "TaxService3", NULL, 30);
}

#ifdef OS_MAC
#include <Events.h>
#endif

NLM_EXTERN Taxon3ReplyPtr Tax3WaitForReply (
  CONN conn
)

{
  AsnIoConnPtr    aicp;
  time_t          currtime, starttime;
  time_t          max = 0;
  EIO_Status      status;
  STimeout        timeout;
  Taxon3ReplyPtr  t3ry = NULL;
#ifdef OS_MAC
  EventRecord     currEvent;
#endif

  if (conn == NULL) return NULL;

#ifdef OS_MAC
  timeout.sec = 0;
  timeout.usec = 0;
#else
  timeout.sec = 300;
  timeout.usec = 0;
#endif

  starttime = GetSecs ();
  while ((status = CONN_Wait (conn, eIO_Read, &timeout)) == eIO_Timeout && max < 300) {
    currtime = GetSecs ();
    max = currtime - starttime;
#ifdef OS_MAC
    WaitNextEvent (0, &currEvent, 0, NULL);
#endif
  }
  if (status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "r" : "rb", conn);
    t3ry = Taxon3ReplyAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }
  CONN_Close (conn);

  return t3ry;
}

/* high-level connection functions */

NLM_EXTERN Taxon3ReplyPtr Tax3SynchronousQuery (
  Taxon3RequestPtr t3rq
)

{
  AsnIoConnPtr    aicp;
  CONN            conn;
  Taxon3ReplyPtr  t3ry;

  if (t3rq == NULL) return NULL;

  conn = Tax3OpenConnection ();

  if (conn == NULL) return NULL;

  aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "w" : "wb", conn);

  Taxon3RequestAsnWrite (t3rq, aicp->aip, NULL);

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  t3ry = Tax3WaitForReply (conn);

  return t3ry;
}

NLM_EXTERN Boolean Tax3AsynchronousQuery (
  Taxon3RequestPtr t3rq,
  QUEUE* queue,
  QueryResultProc resultproc,
  VoidPtr userdata
)

{
  AsnIoConnPtr  aicp;
  CONN          conn;

  if (t3rq == NULL) return FALSE;

  conn = Tax3OpenConnection ();

  if (conn == NULL) return FALSE;

  aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "w" : "wb", conn);

  Taxon3RequestAsnWrite (t3rq, aicp->aip, NULL);

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (queue, conn, resultproc, userdata, TRUE);

  return TRUE;
}

NLM_EXTERN Int4 Tax3CheckQueue (
  QUEUE* queue
)

{
  return QUERY_CheckQueue (queue);
}

NLM_EXTERN Taxon3ReplyPtr Tax3ReadReply (
  CONN conn,
  EIO_Status status
)

{
  AsnIoConnPtr    aicp;
  Taxon3ReplyPtr  t3ry = NULL;

  if (conn != NULL && status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "r" : "rb", conn);
    t3ry = Taxon3ReplyAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }
  return t3ry;
}

NLM_EXTERN Taxon3RequestPtr CreateTaxon3Request (
  Int4 taxid,
  CharPtr name,
  OrgRefPtr orp
)

{
  Taxon3RequestPtr  t2rp;

  t2rp = Taxon3RequestNew ();
  if (t2rp == NULL) return NULL;

  if (StringDoesHaveText (name)) {
    ValNodeCopyStr (&(t2rp->request), 2, name);
  } else if (taxid > 0) {
    ValNodeAddInt (&(t2rp->request), 1, taxid);
  } else if (orp != NULL) {
    orp = AsnIoMemCopy ((Pointer) orp,
                        (AsnReadFunc) OrgRefAsnRead,
                        (AsnWriteFunc) OrgRefAsnWrite);
    ValNodeAddPointer (&(t2rp->request), 3, (Pointer) orp);
  }

  return t2rp;
}

NLM_EXTERN Taxon3RequestPtr CreateMultiTaxon3Request (ValNodePtr org_list)
{
  ValNodePtr vnp;
  Taxon3RequestPtr t3rp;
  OrgRefPtr orp;
  
  t3rp = Taxon3RequestNew ();
  if (t3rp == NULL) return NULL;

  for (vnp = org_list; vnp != NULL; vnp = vnp->next)
  {
    switch (vnp->choice)
    {
      case 1:
        ValNodeAddInt (&(t3rp->request), 1, vnp->data.intvalue);
        break;
      case 2:
        ValNodeCopyStr (&(t3rp->request), 2, vnp->data.ptrvalue);
        break;
      case 3:
        orp = AsnIoMemCopy (vnp->data.ptrvalue,
                        (AsnReadFunc) OrgRefAsnRead,
                        (AsnWriteFunc) OrgRefAsnWrite);
        ValNodeAddPointer (&(t3rp->request), 3, (Pointer) orp);
        break;
    }
  }
  return t3rp;
}


static Boolean HasMisspellingFlag (T3DataPtr t)
{
  T3StatusFlagsPtr status;

  if (t == NULL) return FALSE;
  status = t->status;
  while (status != NULL) {
    if (StringCmp (status->property, "misspelled_name") == 0) {
      return TRUE;
    }
    status = status->next;
  }
  return FALSE;
}


NLM_EXTERN ValNodePtr Taxon3GetOrgRefList (ValNodePtr org_list)
{
  Taxon3RequestPtr t3rq;
  Taxon3ReplyPtr   t3ry;
  T3DataPtr        tdp;
  OrgRefPtr        t3orp = NULL;
  T3ReplyPtr       trp;
  T3ErrorPtr       tep;
  ValNodePtr       response_list = NULL, next_org_list, last_org;
  Int4             request_num, max_requests = 2000;

  while (org_list != NULL) {
    /* we need to break large org_lists into manageable chunks */
    next_org_list = org_list->next;
    last_org = org_list; 
    request_num = 1;
    while (next_org_list != NULL && request_num < max_requests) {
      last_org = next_org_list;
      next_org_list = next_org_list->next;
      request_num++;
    }
    if (last_org != NULL) {
      last_org->next = NULL;
    }
      
    /* now create the request */
  
    t3rq = CreateMultiTaxon3Request (org_list);
    if (t3rq == NULL) return NULL;
    t3ry = Tax3SynchronousQuery (t3rq);
    Taxon3RequestFree (t3rq);
    if (t3ry != NULL) {
      for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
        switch (trp->choice) {
          case T3Reply_error :
            tep = (T3ErrorPtr) trp->data.ptrvalue;
            if (tep != NULL) {
              ErrPostEx (SEV_ERROR, 0, 0, tep->message);
            }
            ValNodeAddPointer (&response_list, 3, NULL);
            break;
          case T3Reply_data :
            tdp = (T3DataPtr) trp->data.ptrvalue;
            if (tdp != NULL) {
              t3orp = (OrgRefPtr)(tdp->org);
              if (HasMisspellingFlag (tdp)) {
                ValNodeAddPointer (&response_list, 4, (Pointer) t3orp);
              } else {
                ValNodeAddPointer (&response_list, 3, (Pointer) t3orp);
              }
              tdp->org = NULL;
            }
            break;
          default :
            break;
        }
      }
      Taxon3ReplyFree (t3ry);
    }
    
    if (last_org != NULL) {
        last_org->next = next_org_list;
    }
    org_list = next_org_list;
  }
  return response_list;
}

NLM_EXTERN OrgRefPtr Taxon3GetOrg (OrgRefPtr orp)

{
  Taxon3RequestPtr t3rq;
  Taxon3ReplyPtr   t3ry;
  T3DataPtr        tdp;
  OrgRefPtr        t3orp = NULL;
  T3ReplyPtr        trp;
  T3ErrorPtr        tep;
	
  if (orp == NULL) return NULL;
  
  t3rq = CreateTaxon3Request (0, NULL, orp);
  if (t3rq == NULL) return NULL;
  t3ry = Tax3SynchronousQuery (t3rq);
  Taxon3RequestFree (t3rq);
  if (t3ry != NULL) {
    for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
      switch (trp->choice) {
        case T3Reply_error :
          tep = (T3ErrorPtr) trp->data.ptrvalue;
          if (tep != NULL) {
            ErrPostEx (SEV_ERROR, 0, 0, tep->message);
          }
          break;
        case T3Reply_data :
          tdp = (T3DataPtr) trp->data.ptrvalue;
          if (tdp != NULL) {
            t3orp = (OrgRefPtr)(tdp->org);
            tdp->org = NULL;
          }
          break;
        default :
          break;
      }
    }
    Taxon3ReplyFree (t3ry);
  }
  
  return t3orp;
}

static Boolean DoOrgIdsMatch(BioSourcePtr b1, BioSourcePtr b2)
{
  DbtagPtr d1 = NULL, d2 = NULL;
  ValNodePtr vnp;
	
  if (b1 == NULL || b2 == NULL) 
  {
    return FALSE;
  }
  if (b1->org ==  NULL || b2->org == NULL) 
  {
    return FALSE;
  }
  for (vnp = b1->org->db; vnp; vnp = vnp->next) 
  {
    d1 = (DbtagPtr) vnp->data.ptrvalue;
    if (StringCmp(d1->db, "taxon") == 0) 
    {
      break;
    }
  }
  for (vnp = b2->org->db; vnp; vnp = vnp->next) 
  {
    d2 = (DbtagPtr) vnp->data.ptrvalue;
	if (StringCmp(d2->db, "taxon") == 0) 
	{
      break;
	}
  }
  if (d1 && d2) 
  {
	if (d1->tag->id == d2->tag->id) 
	{
      return TRUE;
	}
  }
  else if (StringICmp(b1->org->taxname, b2->org->taxname) == 0) 
  {
	return TRUE;
  }
  return FALSE;
}

static BioSourcePtr Tax3BioSourceMerge(BioSourcePtr host, BioSourcePtr guest)
{
  SubSourcePtr ssp, sp, last_ssp;
  OrgModPtr omp, homp, last_omp;
  OrgNamePtr	onp;
	
  if (host == NULL && guest == NULL) 
  {
    return NULL;
  }
  if (host == NULL && guest != NULL) 
  {
	host = AsnIoMemCopy(guest, (AsnReadFunc) BioSourceAsnRead, 
		   						(AsnWriteFunc) BioSourceAsnWrite);
	return host;
  }
  if (host != NULL && guest == NULL) 
  {
    return host;
  }
  if (host->genome == 0 && guest->genome != 0) 
  {
    host->genome = guest->genome;
  }
  if (host->origin == 0 && guest->origin != 0) 
  {
    host->origin = guest->origin;
  }
  last_ssp = host->subtype;
  while (last_ssp != NULL && last_ssp->next != NULL)
  {
  	last_ssp = last_ssp->next;
  }
  for (ssp = guest->subtype; ssp; ssp = ssp->next) 
  {
    sp = AsnIoMemCopy(ssp, (AsnReadFunc) SubSourceAsnRead, 
		   						(AsnWriteFunc) SubSourceAsnWrite);
    if (last_ssp == NULL)
    {
      host->subtype = sp;
    }
    else
    {
      last_ssp->next = sp;
      last_ssp = sp;
    }
  }
  if (guest->org->orgname) 
  {
   	if ((onp = host->org->orgname)	== NULL) 
   	{
   	  onp = OrgNameNew();
   	  host->org->orgname = onp;
    }	
    last_omp = onp->mod;		
    while (last_omp != NULL && last_omp->next != NULL)
    {
      last_omp = last_omp->next;
    }
    for (omp = guest->org->orgname->mod; omp; omp = omp->next) 
    {
      homp = AsnIoMemCopy(omp, (AsnReadFunc) OrgModAsnRead, 
		   						(AsnWriteFunc) OrgModAsnWrite);
      if (last_omp == NULL)
      {
      	onp->mod = homp;
      }
      else
      {
      	last_omp->next = homp;
      	last_omp = homp;
      }
    }
  }
  return host;
}


/**************************************************************************
*	Compare BioSources in one bioseq->descr using Taxonomy to find
*	their join parent
*	merge if organisms are the same or create a feature if different
*
**************************************************************************/
NLM_EXTERN void Tax3MergeSourceDescr (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr    bsp = NULL;
	ValNodePtr   vnp, newlist;
	SeqFeatPtr   sfp;
	BioSourcePtr first_biop = NULL;
	BioSourcePtr other_biop;
	BioSourcePtr tmp_biop;
	ObjValNodePtr ovp;

	if (!IS_Bioseq(sep)) {
		return;
	}
	newlist = (ValNodePtr) data;
	bsp = (BioseqPtr) sep->data.ptrvalue;
	if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const) 
			&& (bsp->repr != Seq_repr_delta))
		return;

	if (! ISA_na(bsp->mol))
		return;
	
	/* add the descriptors in newlist to the end of the list in bsp->descr*/
	if (bsp->descr == NULL)
	{
	  bsp->descr = newlist;
	}
	else
	{
	  for (vnp = bsp->descr; vnp->next != NULL; vnp = vnp->next)
	  {	
	  }
	  vnp->next = newlist;
	}
	
	/* now find the first source descriptor in bsp->descr that has an org*/
    /* note - we can't use SeqMgrGetNextDescriptor here because we have just
     * added to the descriptors, so they are not indexed. */
	for (vnp = bsp->descr; vnp != NULL; vnp = vnp->next)
	{
	  if (vnp->choice != Seq_descr_source) continue;
	  if (vnp->data.ptrvalue == NULL)
	  {
	  	ErrPostStr(SEV_WARNING, 0, 0, "Source descriptor missing data");
	  	if (vnp->extended)
	  	{
	  	  ovp = (ObjValNodePtr) vnp;
	  	  ovp->idx.deleteme = TRUE;
	  	}
	  }
	  if (first_biop == NULL)
	  {
	  	first_biop = vnp->data.ptrvalue;
	  }
	  else
	  {
		other_biop = vnp->data.ptrvalue;
		/* detach biosource pointer from descr, so that it will not be freed
		 * when the descriptor is deleted.
		 */
		vnp->data.ptrvalue = NULL;
        if (vnp->extended)
        {
          ovp = (ObjValNodePtr) vnp;
	  	  ovp->idx.deleteme = TRUE;
        }
        if (DoOrgIdsMatch(first_biop, other_biop)) 
		{
		  /* merge the two sources */
		  tmp_biop = Tax3BioSourceMerge(first_biop, other_biop);
		  if (tmp_biop == NULL)
		  {
		  	ErrPostStr (SEV_WARNING, 0, 0, "Failed to merge biosources");
		  }
		  else
		  {
		  	first_biop = tmp_biop;
		  }
		  other_biop = BioSourceFree (other_biop);
		} else {
		  /* create a source feature */
		  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_BIOSRC, NULL);
		  if (sfp != NULL)
		  {
            sfp->data.value.ptrvalue = other_biop;
		  }
        }
	  }
	}
	return;
}

static Int4 GetTaxIdFromOrgRef (OrgRefPtr orp)
{
  Int4       tax_id = -1;
  ValNodePtr vnp;
  DbtagPtr   d;

  if (orp != NULL)
  {
    for (vnp = orp->db; vnp != NULL; vnp = vnp->next) 
    {
      d = (DbtagPtr) vnp->data.ptrvalue;
      if (StringCmp(d->db, "taxon") == 0) 
      {
        tax_id = d->tag->id;
        break;
      }
    }
  }
  return tax_id;
}

NLM_EXTERN Int4 Taxon3GetTaxIdByOrgRef (OrgRefPtr orp)
{
  OrgRefPtr  orp_repl;
  Int4       tax_id = -1;
  
  if (orp == NULL) return -1;
  
  orp_repl = Taxon3GetOrg (orp);
  tax_id = GetTaxIdFromOrgRef (orp_repl);
  OrgRefFree (orp_repl);
  
  return tax_id;
}

NLM_EXTERN OrgRefPtr Taxon3GetOrgRefByName (CharPtr orgname)
{
  OrgRefPtr request, org;
  
  request = OrgRefNew ();
  if (request == NULL) return NULL;
  request->taxname = orgname;
  org = Taxon3GetOrg (request);
  request->taxname = NULL;
  OrgRefFree (request);
  return org;
}

NLM_EXTERN Int4 Taxon3GetTaxIdByName (CharPtr orgname)
{
  OrgRefPtr orp;
  Int4      tax_id;
  
  orp = Taxon3GetOrgRefByName (orgname);
  tax_id = GetTaxIdFromOrgRef (orp);

  OrgRefFree(orp);
  return tax_id;
}

static void AddBioSourceToList (BioSourcePtr biop, Pointer userdata)
{
  ValNodePtr PNTR list;
  
  if (biop == NULL || userdata == NULL) return;
  list = (ValNodePtr PNTR) userdata;
  ValNodeAddPointer (list, 4, (Pointer) biop);
}

NLM_EXTERN void Taxon3ReplaceOrgInSeqEntry (SeqEntryPtr sep, Boolean keep_syn)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   biop_vnp, response_vnp;
  BioSourcePtr biop;
  OrgRefPtr    swap_org, response_org;
  
  VisitBioSourcesInSep (sep, &biop_list, AddBioSourceToList);

  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {
    biop = (BioSourcePtr) biop_vnp->data.ptrvalue;
    ValNodeAddPointer (&request_list, 3, biop->org);
  }
  response_list = Taxon3GetOrgRefList (request_list);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
    return;
  }

  for (biop_vnp = biop_list, response_vnp = response_list;
       biop_vnp != NULL && response_vnp != NULL;
       biop_vnp = biop_vnp->next, response_vnp = response_vnp->next)
  {
    biop = (BioSourcePtr) biop_vnp->data.ptrvalue;
    swap_org = biop->org;
    response_org = response_vnp->data.ptrvalue;
    if (response_org != NULL)
    {
      biop->org = response_org;
      response_vnp->data.ptrvalue = NULL;
      OrgRefFree (swap_org);
      if (! keep_syn)
      {
        biop->org->syn = ValNodeFreeData(biop->org->syn);
      }
    }
  }
  ValNodeFree (request_list);
  ValNodeFree (response_list);
  ValNodeFree (biop_list);   
}


static void GetBioSourceFeaturesForCheck (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR list = (ValNodePtr PNTR) userdata;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || list == NULL
      || sfp->data.value.ptrvalue == NULL) {
    return;
  }
  ValNodeAddPointer (list, OBJ_SEQFEAT, sfp);
}


static void GetBioSourceDescriptorsForCheck (SeqDescrPtr sdp, Pointer userdata)
{
  ValNodePtr PNTR list = (ValNodePtr PNTR) userdata;
  if (sdp == NULL || sdp->choice != Seq_descr_source || list == NULL
      || sdp->data.ptrvalue == NULL) {
    return;
  }
  ValNodeAddPointer (list, OBJ_SEQDESC, sdp);
}


static DbtagPtr GetTaxonXref (OrgRefPtr org)
{
  ValNodePtr vnp;
  DbtagPtr   dbt = NULL;
  
  if (org == NULL) return NULL;
  vnp = org->db;
  while (vnp != NULL && dbt == NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && StringICmp ((CharPtr) dbt->db, "taxon") != 0) {
      dbt = NULL;
    }
    vnp = vnp->next;
  }
  return dbt;
}
  
static Boolean DoTaxonIdsMatch (OrgRefPtr org1, OrgRefPtr org2)
{
  DbtagPtr   dbt1 = NULL, dbt2 = NULL;
  
  if (org1 == NULL || org2 == NULL) return FALSE;
  
  dbt1 = GetTaxonXref (org1);
  if (dbt1 == NULL) return FALSE;
  dbt2 = GetTaxonXref (org2);
  if (dbt2 == NULL) return FALSE;
  
  return DbtagMatch(dbt1, dbt2);
}


NLM_EXTERN void Taxon3CheckOrgInSeqEntry (SeqEntryPtr sep, ValNodePtr PNTR not_found, ValNodePtr PNTR bad_match)
{
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   biop_vnp, response_vnp;
  BioSourcePtr biop;
  OrgRefPtr    orig_org, response_org;
  ValNodePtr   item_list = NULL;
  SeqFeatPtr   sfp;
  SeqDescrPtr  sdp;
  
  VisitFeaturesInSep (sep, &item_list, GetBioSourceFeaturesForCheck);
  VisitDescriptorsInSep (sep, &item_list, GetBioSourceDescriptorsForCheck);
  
  for (biop_vnp = item_list; biop_vnp != NULL; biop_vnp = biop_vnp->next) {
    biop = NULL;
    if (biop_vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) biop_vnp->data.ptrvalue;  
      if (sfp != NULL) {  
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;      
      }
    } else if (biop_vnp->choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) biop_vnp->data.ptrvalue;
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      }
    }
    if (biop != NULL) {
      ValNodeAddPointer (&request_list, 3, biop->org);
    }
  }

  response_list = Taxon3GetOrgRefList (request_list);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
    ValNodeFree (request_list);
    ValNodeFree (item_list);
    return;
  }

  for (biop_vnp = item_list, response_vnp = response_list;
       biop_vnp != NULL && response_vnp != NULL;
       biop_vnp = biop_vnp->next, response_vnp = response_vnp->next)
  {
    response_org = response_vnp->data.ptrvalue;  
    biop = NULL;
    orig_org = NULL;
    if (biop_vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) biop_vnp->data.ptrvalue;    
      if (sfp != NULL) {  
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      }
    } else if (biop_vnp->choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) biop_vnp->data.ptrvalue;
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      }
    }
    if (biop == NULL) {
      Message (MSG_POST, "Error collecting data");
      ValNodeFree (request_list);
      ValNodeFree (item_list);
      return;
    } else {
      orig_org = biop->org;
      if (orig_org != NULL) {
        if (response_org == NULL) {
          ValNodeAddPointer (not_found, biop_vnp->choice, biop_vnp->data.ptrvalue);          
        } else if (StringCmp (orig_org->taxname, response_org->taxname) != 0) {
          ValNodeAddPointer (bad_match, biop_vnp->choice, biop_vnp->data.ptrvalue);
        } else if (!DoTaxonIdsMatch(orig_org, response_org)) {
          ValNodeAddPointer (bad_match, biop_vnp->choice, biop_vnp->data.ptrvalue);
        }        
      }
    }
    OrgRefFree (response_org);
  }
  ValNodeFree (request_list);
  ValNodeFree (response_list);
  ValNodeFree (item_list);   
}


NLM_EXTERN void CheckTaxNamesAgainstTaxDatabase (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp;
  SeqEntryPtr sep;
  SeqEntryPtr orig_scope;
  ValNodePtr  not_found = NULL, bad_match = NULL;
  CharPtr     bad_match_fmt = "%d tax names do not match taxonomy lookup.";
  CharPtr     no_match_fmt = "%d organisms are not found in taxonomy lookup.";
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  
  orig_scope = SeqEntryGetScope ();
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    SeqEntrySetScope (sep);
    Taxon3CheckOrgInSeqEntry (sep, &not_found, &bad_match);
  }
  SeqEntrySetScope (orig_scope);
  if (not_found != NULL) {
    dip = NewClickableItem (DISC_NO_TAXLOOKUP, no_match_fmt, not_found);
    dip->subcategories = NULL;
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
  if (bad_match != NULL) {
    dip = NewClickableItem (DISC_BAD_TAXLOOKUP, bad_match_fmt, bad_match);
    dip->subcategories = NULL;
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static ValNodePtr FreeOrgRefValNodeList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  OrgRefPtr  org;

  while (vnp != NULL)
  { 
    vnp_next = vnp->next;
    vnp->next = NULL;
    org = (OrgRefPtr) vnp->data.ptrvalue;
    vnp->data.ptrvalue = OrgRefFree (org);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static Boolean EndsWithSp (CharPtr str)
{
  Int4 len;

  if (StringHasNoText (str)) return FALSE;
  len = StringLen (str);
  if (len < 4) return FALSE;
  if (StringCmp (str + len - 4, " sp.") == 0) return TRUE;
  return FALSE;
}


static CharPtr RemoveSp (CharPtr orig)
{
  CharPtr cpy = NULL;
  Int4    len;

  len = StringLen (orig);
  if (len >= 4 && StringCmp (orig + len - 4, " sp.") == 0) {
    cpy = (CharPtr) MemNew (sizeof (Char) * len - 3);
    StringNCpy (cpy, orig, len - 4);
    cpy[len - 4] = 0;
  }
  return cpy;
}

  
static void AddRequestOrgForString (CharPtr str, CharPtr host, ValNodePtr PNTR request_list, ValNodePtr PNTR req_host_list)
{
  OrgRefPtr    request_org;
  CharPtr      cp, cpy;
  Int4         word1_len, space_len, word2_len;

  if (StringHasNoText (str) || host == NULL || request_list == NULL || req_host_list == NULL)
  {
    return;
  }

  /* if ends with " sp.", remove " sp." */
  cpy = RemoveSp (host);
  if (cpy != NULL) {
    request_org = OrgRefNew();
    request_org->taxname = cpy;
    ValNodeAddPointer (request_list, 3, request_org);
    ValNodeAddPointer (req_host_list, 0, host);
  } else {
    request_org = OrgRefNew();
    request_org->taxname = StringSave (str);
    ValNodeAddPointer (request_list, 3, request_org);
    ValNodeAddPointer (req_host_list, 0, host);

     
    /* if more than one word, try chopping off last to see if abbreviated name looks up */
    word1_len = StringCSpn (str, " ");
    if (word1_len == 0) return;
    space_len = StringSpn (str + word1_len, " ");
    if (space_len == 0) return;
    word2_len = StringCSpn (str + word1_len + space_len, " ");
    if (word2_len == 0) return;
    if (isspace (*(str + word1_len + space_len + word2_len)))
    {
      cpy = StringSave (str);    
      cp = StringRChr (cpy, ' ');
      if (cp != NULL)
      {
        *cp = 0;
        AddRequestOrgForString (cpy, host, request_list, req_host_list);
      }
      cpy = MemFree (cpy);
    }
  }
}

typedef struct specifichostcheck {
  CharPtr      spec_host;
  ValNodePtr   request_list;  /* ValNodeList of orgs */
  ValNodePtr   response_list; /* ValNodeList of orgs */
  ValNodePtr   biop_list;     /* ValNodeList of sources with this spec_host value */
} SpecificHostCheckData, PNTR SpecificHostCheckPtr;


static ValNodePtr SpecificHostCheckListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  SpecificHostCheckPtr p;

  while (vnp != NULL)
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    p = (SpecificHostCheckPtr) vnp->data.ptrvalue;
    if (p != NULL)
    {
      p->request_list = FreeOrgRefValNodeList (p->request_list);
      p->response_list = FreeOrgRefValNodeList (p->response_list);
      p->spec_host = MemFree (p->spec_host);
      p->biop_list = ValNodeFree (p->biop_list);
    }
    vnp = ValNodeFreeData (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static ValNodePtr SortSpecificHostOrgs (ValNodePtr host_list, ValNodePtr request_list, ValNodePtr response_list)
{
  ValNodePtr           check_list = NULL;
  SpecificHostCheckPtr p = NULL;
  CharPtr              host, prev_host = NULL;

  while (host_list != NULL
         && request_list != NULL
         && response_list != NULL)
  {
    host = (CharPtr) host_list->data.ptrvalue;
    if (StringCmp (host, prev_host) != 0)
    {
      p = (SpecificHostCheckPtr) MemNew (sizeof (SpecificHostCheckData));
      p->spec_host = StringSave (host);
      ValNodeAddPointer (&check_list, 0, p);
      prev_host = host;
    }
    ValNodeAddPointer (&(p->request_list), request_list->choice, request_list->data.ptrvalue);
    ValNodeAddPointer (&(p->response_list), response_list->choice, response_list->data.ptrvalue);
    request_list->data.ptrvalue = NULL;
    response_list->data.ptrvalue = NULL;
    host_list = host_list->next;
    request_list = request_list->next;
    response_list = response_list->next;
  }
  return check_list;        
}


static Boolean StringAlreadyInValNodeList (CharPtr str, ValNodePtr list) 
{
  if (StringHasNoText (str))
  {
    return TRUE;
  }
  
  while (list != NULL)
  {
    if (StringCmp (str, list->data.ptrvalue) == 0)
    {
      return TRUE;
    }
    list = list->next;
  }
  return FALSE;
}


static BioSourcePtr GetBioSourceFromValNode (ValNodePtr vnp)
{
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioSourcePtr biop = NULL;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;

  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  } 
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }
  return biop;
}

static void AddBioSourcesToSpecificHostChecklist (ValNodePtr biop_list, ValNodePtr check_list)
{
  ValNodePtr biop_vnp, last_vnp = NULL, stop_search;
  BioSourcePtr biop;
  OrgModPtr    mod;
  SpecificHostCheckPtr p;

  if (biop_list == NULL || check_list == NULL) return;

  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {

    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL) continue;

    if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) continue;
    mod = biop->org->orgname->mod;
    while (mod != NULL)
    {
      if (mod->subtype == ORGMOD_nat_host
          && !StringHasNoText (mod->subname))
      {
        if (last_vnp == NULL)
        {
          last_vnp = check_list;
          stop_search = NULL;
        }
        else
        {
          stop_search = last_vnp;
        }
        p = NULL;
        while (last_vnp != NULL 
               && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
               && StringCmp (p->spec_host, mod->subname) != 0)
        {
          p = NULL;
          last_vnp = last_vnp->next;
        }
        if (p == NULL && stop_search != NULL)
        {
          last_vnp = check_list;
          while (last_vnp != stop_search 
                 && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
                 && StringCmp (p->spec_host, mod->subname) != 0)
          {
            p = NULL;
            last_vnp = last_vnp->next;
          }
        }

        if (p != NULL)
        {
          ValNodeAddPointer (&(p->biop_list), biop_vnp->choice, biop_vnp->data.ptrvalue);
        }
      }
      mod = mod->next;
    }
  }
}


NLM_EXTERN Boolean IsOrgModSpecificHostToBeChecked (OrgModPtr mod, Boolean for_validator, Boolean check_single_word)
{
  CharPtr   cp;
  Boolean   any_upper = FALSE, any_space = FALSE;

  if (mod != NULL && mod->subtype == ORGMOD_nat_host
      && !StringHasNoText (mod->subname))
  {
    cp = mod->subname;
    while (*cp != 0 && (!any_upper || !any_space))
    {
      if (isupper (*cp))
      {
        any_upper = TRUE;
      }
      else if (*cp == ' ')
      {
        any_space = TRUE;
      }
      else if (ispunct (*cp))
      {
        return TRUE;
      }
      cp++;
    }
    if (any_upper && (any_space || check_single_word))
    {
      return TRUE;
    }
  }
  return FALSE;
}


static Boolean HasSpecificHostToBeChecked (BioSourcePtr biop, Boolean for_validator, Boolean check_single_word)
{
  OrgModPtr mod;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) return FALSE;
  mod = biop->org->orgname->mod;
  while (mod != NULL)
  {
    if (IsOrgModSpecificHostToBeChecked (mod, for_validator, check_single_word))
    {
      return TRUE;
    }
    mod = mod->next;
  }
  return FALSE;
}


static Boolean MatchesSynonym (CharPtr txt, OrgRefPtr response_org)
{
  ValNodePtr syn;
  Boolean    rval = FALSE;
  if (StringHasNoText (txt) || response_org == NULL) return FALSE;

  for (syn = response_org->syn; syn != NULL && !rval; syn = syn->next)
  {
    if (StringCmp (txt, syn->data.ptrvalue) == 0)
    {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean MatchesGenBankSynonym (CharPtr txt, OrgRefPtr response_org)
{
  OrgModPtr mod;
  Boolean   rval = FALSE;

  if (StringHasNoText (txt) || response_org == NULL || response_org->orgname == NULL) return FALSE;
  mod = response_org->orgname->mod;
  while (mod != NULL) 
  {
    if ((mod->subtype == ORGMOD_gb_synonym || mod->subtype == ORGMOD_old_name) && StringCmp (txt, mod->subname) == 0)
    {
      rval = TRUE;
    }
    mod = mod->next;
  }
  return rval;
}


static Boolean MatchesCommonName (CharPtr txt, CharPtr common_name)
{
  CharPtr cp;
  Int4    len;

  if (StringHasNoText (txt) || StringHasNoText (common_name))
  {
    return FALSE;
  }
  else if (StringICmp (txt, common_name) == 0)
  {
    return TRUE;
  }
  else
  {
    cp = StringISearch (txt, common_name);
    len = StringLen (common_name);
    if (cp != NULL
        && (cp == txt || isspace (*(cp - 1)))
        && (*(cp + len) == 0 || isspace (*(cp + len))))
    {
      return TRUE;
    }
    cp = StringISearch (common_name, txt);
    len = StringLen (txt);
    if (cp != NULL
        && (cp == common_name || isspace (*(cp - 1)))
        && (*(cp + len) == 0 || isspace (*(cp + len))))
    {
      return TRUE;
    }
  }

  return FALSE;
}

typedef struct spechostgather {
  ValNodePtr list;
  Boolean    for_validator;
  Boolean    check_single_word;
} SpecHostGatherData, PNTR SpecHostGatherPtr;

static void AddSpecificHostBioSourceFeatToList (SeqFeatPtr sfp, Pointer userdata)
{
  SpecHostGatherPtr p;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;

  p = (SpecHostGatherPtr) userdata;
  if (HasSpecificHostToBeChecked (sfp->data.value.ptrvalue, p->for_validator, p->check_single_word))
  {
    ValNodeAddPointer (&(p->list), OBJ_SEQFEAT, sfp);
  }
}


static void AddSpecificHostBioSourceDescToList (SeqDescrPtr sdp, Pointer userdata)
{
  SpecHostGatherPtr p;

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;

  p = (SpecHostGatherPtr) userdata;
  if (HasSpecificHostToBeChecked (sdp->data.ptrvalue, p->for_validator, p->check_single_word))
  {
    ValNodeAddPointer (&(p->list), OBJ_SEQDESC, sdp);
  }
}


static ValNodePtr GetSpecificHostBioSourceList (SeqEntryPtr sep, Boolean for_validator, Boolean check_single_word)
{
  SpecHostGatherData   d;
  
  d.for_validator = for_validator;
  d.check_single_word = check_single_word;
  d.list = NULL;
  VisitFeaturesInSep (sep, &d, AddSpecificHostBioSourceFeatToList);
  VisitDescriptorsInSep (sep, &d, AddSpecificHostBioSourceDescToList);
  return d.list;
}

static ValNodePtr GetListOfUniqueSpecificHostValues (ValNodePtr biop_list)
{
  ValNodePtr   biop_vnp;
  BioSourcePtr biop;
  OrgModPtr    mod;
  ValNodePtr   spec_host_list = NULL;
  
  /* get a list of unique specific_host values */
  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {
    if (biop_vnp->data.ptrvalue == NULL) continue;
    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) continue;
    mod = biop->org->orgname->mod;
    while (mod != NULL)
    {
      if (mod->subtype == ORGMOD_nat_host
          && !StringHasNoText (mod->subname)
          && !StringAlreadyInValNodeList (mod->subname, spec_host_list))
      {
        ValNodeAddPointer (&spec_host_list, 0, mod->subname);
      }
      mod = mod->next;
    }
  }
  return spec_host_list;
}


static void 
FormatSpecificHostRequests 
(ValNodePtr spec_host_list,
 ValNodePtr PNTR request_list,
 ValNodePtr PNTR req_host_list)
{
  ValNodePtr vnp;
  CharPtr    orig, cp, str, cp2 = NULL;
  
  /* now format requests for unique specific_host values */
  for (vnp = spec_host_list; vnp != NULL; vnp = vnp->next)
  {
    orig = (CharPtr) vnp->data.ptrvalue;
    /* if we have a value in parentheses, submit it separately */
    cp = StringChr (orig, '(');
    if (cp != NULL)
    {
      cp2 = StringChr (cp, ')');
    }
    if (cp != NULL && cp2 != NULL 
        && ((cp > orig && orig[StringLen (orig) - 1] == ')') /* ends with paren */
            || (cp == orig))) /* starts with paren */
    {
      if (cp > orig && orig[StringLen (orig) - 1] == ')')
      {
        str = StringSave (orig);
        /* remove trailing parenthesis */
        str [StringLen(str) - 1] = 0;

        cp = str + (cp - orig);

        /* remove opening parenthesis */
        *cp = 0;
        cp++;
      }
      else
      {
        str = StringSave (orig);
        /* remove leading parenthesis */
        str[0] = ' ';
        cp = str + (cp2 - orig);
        /* remove trailing parenthesis */
        *cp = 0; 
        cp++;
      }
      TrimSpacesAroundString (cp);
      TrimSpacesAroundString (str);
      AddRequestOrgForString (str, orig, request_list, req_host_list);
      AddRequestOrgForString (cp, orig, request_list, req_host_list);
    }
    else
    {
      AddRequestOrgForString (orig, orig, request_list, req_host_list);
    }
  }
}


static Boolean IsShorterVersionOfEarlierMatch (ValNodePtr request_vnp, ValNodePtr request_list, ValNodePtr response_list)
{
  OrgRefPtr this_org, earlier_org, earlier_response;
  Int4      len;
  Boolean   found = FALSE;

  if (request_vnp == NULL || request_list == NULL) return FALSE;
  this_org = (OrgRefPtr) request_vnp->data.ptrvalue;
  if (this_org == NULL || StringHasNoText (this_org->taxname)) return FALSE;
  len = StringLen (this_org->taxname);
  while (request_list != NULL && request_list != request_vnp && !found && response_list != NULL)
  {
    earlier_org = (OrgRefPtr) request_list->data.ptrvalue;
    earlier_response = (OrgRefPtr) response_list->data.ptrvalue;
    if (earlier_org != NULL && earlier_response != NULL
        && StringNCmp (earlier_org->taxname, this_org->taxname, len) == 0
        && (MatchesCommonName (earlier_org->taxname, earlier_response->common)
            || MatchesSynonym (earlier_org->taxname, earlier_response)
            || MatchesGenBankSynonym (earlier_org->taxname, earlier_response)))
    { 
      found = TRUE;
    }

    request_list = request_list->next;
    response_list = response_list->next;
  }
  return found;
}


static Boolean MatchWithSp (CharPtr request, CharPtr response)
{
  CharPtr cpy = NULL;
  Boolean rval = FALSE;

  if (StringCmp (request, response) == 0) {
    rval = TRUE;
  } else {
    cpy = RemoveSp (request);
    if (cpy != NULL && StringCmp (cpy, response) == 0) {
      rval = TRUE;
    }
    cpy = MemFree (cpy);
  }
  return rval;
}


static Boolean FindMatchInResponseOrgList (CharPtr request_str, ValNodePtr response_list)
{
  ValNodePtr vnp;
  OrgRefPtr  org;
  Boolean    rval = FALSE;

  if (StringHasNoText (request_str) || response_list == NULL) {
    return FALSE;
  }

  for (vnp = response_list; vnp != NULL && !rval; vnp = vnp->next) {
    org = (OrgRefPtr) vnp->data.ptrvalue;
    if (org != NULL) {
      if (MatchesCommonName (request_str, org->common)
          || MatchesSynonym (request_str, org)) {
        rval = TRUE;
      }
    }
  }
  return rval;
}


/* Want to check that specific host names are valid */
NLM_EXTERN ValNodePtr Taxon3CheckSpecificHostInSeqEntry (SeqEntryPtr sep, Boolean for_validator, Boolean check_single_word)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   req_host_list = NULL, spec_host_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   response_vnp, request_vnp;
  ValNodePtr   check_list, check_vnp;
  OrgRefPtr    request_org, response_org;
  SpecificHostCheckPtr p;
  Boolean              has_taxname, has_bad, all_synonyms;
  ValNodePtr           bad_biop_list = NULL;
  ErrSev               level;
    
  biop_list = GetSpecificHostBioSourceList (sep, for_validator, check_single_word);

  /* get a list of unique specific_host values */
  spec_host_list = GetListOfUniqueSpecificHostValues (biop_list);

  /* now format requests for unique specific_host values */
  FormatSpecificHostRequests (spec_host_list, &request_list, &req_host_list);

  spec_host_list = ValNodeFree (spec_host_list);

  level = ErrSetMessageLevel (SEV_MAX);
  response_list = Taxon3GetOrgRefList (request_list);
  ErrSetMessageLevel (level);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
  }
  else
  {
    /* resort requests so that we can check all responses for the same BioSource together */
    check_list = SortSpecificHostOrgs (req_host_list, request_list, response_list);
    AddBioSourcesToSpecificHostChecklist (biop_list, check_list);  

    /* now look at responses */
    check_vnp = check_list;
    while (check_vnp != NULL)
    {
      p = (SpecificHostCheckPtr) check_vnp->data.ptrvalue;
      if (p != NULL)
      {
        has_taxname = FALSE; 
        has_bad = FALSE;
        request_vnp = p->request_list;
        response_vnp = p->response_list;
        all_synonyms = TRUE;
        while (!has_taxname && request_vnp != NULL && response_vnp != NULL)
        {
          request_org = (OrgRefPtr) request_vnp->data.ptrvalue;
          response_org = (OrgRefPtr) response_vnp->data.ptrvalue;
          if (response_vnp->choice == 4) 
          {
            has_bad = TRUE;
          }  
          else if (response_org == NULL)
          {
            if (!IsShorterVersionOfEarlierMatch(request_vnp, p->request_list, p->response_list)
                && !FindMatchInResponseOrgList (request_org->taxname, p->response_list))
            {
              has_bad = TRUE;
            }
          }
          else if (MatchWithSp (request_org->taxname, response_org->taxname))
          {
            has_taxname = TRUE;
          }          
          else if (MatchesCommonName (request_org->taxname, response_org->common)
                   || MatchesSynonym (request_org->taxname, response_org)
                   || MatchesGenBankSynonym (request_org->taxname, response_org))
          {
            /* it's a synonym */
          }
          else
          {
            has_bad = TRUE;
            all_synonyms = FALSE;
          }
          request_vnp = request_vnp->next;
          response_vnp = response_vnp->next;
        }     
        if (!has_taxname && has_bad)
        {
          /* add to the list of bad */
          ValNodeLink (&bad_biop_list, p->biop_list);
          p->biop_list = NULL;
        }
      }
      check_vnp = check_vnp->next;
    }
    check_list = SpecificHostCheckListFree (check_list);
  }

  biop_list = ValNodeFree (biop_list);
  request_list = FreeOrgRefValNodeList (request_list);
  response_list = FreeOrgRefValNodeList (response_list);
  return bad_biop_list;
}


static Boolean FixSpecificHostForOneBioSource (BioSourcePtr biop, CharPtr new_spec_host, CharPtr old_spec_host)
{
  OrgModPtr mod;
  Boolean   rval = FALSE;

  if (biop == NULL 
      || biop->org == NULL
      || biop->org->orgname == NULL
      || StringHasNoText (new_spec_host) || StringHasNoText (old_spec_host))
  {
    return rval;
  }

  mod = biop->org->orgname->mod;
  while (mod != NULL)
  {
    if (mod->subtype == ORGMOD_nat_host && StringCmp (old_spec_host, mod->subname) == 0)
    {
      mod->subname = MemFree (mod->subname);
      mod->subname = StringSave (new_spec_host);
      rval = TRUE;
    }
    mod = mod->next;
  }
  return rval;
}

static SpecificHostFixPtr SpecificHostFixNew (ValNodePtr feat_or_desc, CharPtr bad_host, CharPtr old_taxname, CharPtr new_taxname)
{
  SpecificHostFixPtr s;

  s = (SpecificHostFixPtr) MemNew (sizeof (SpecificHostFixData));
  if (feat_or_desc != NULL) 
  {
    s->feat_or_desc = ValNodeNew(NULL);
    s->feat_or_desc->choice = feat_or_desc->choice;
    s->feat_or_desc->data.ptrvalue = feat_or_desc->data.ptrvalue;
  }
  s->bad_specific_host = StringSave (bad_host);
  s->old_taxname = StringSave (old_taxname);
  s->new_taxname = StringSave (new_taxname);
  return s;
}


static SpecificHostFixPtr SpecificHostFixFree (SpecificHostFixPtr s)
{
  if (s != NULL)
  {
    s->feat_or_desc = ValNodeFree (s->feat_or_desc);
    s->bad_specific_host = MemFree (s->bad_specific_host);
    s->old_taxname = MemFree (s->old_taxname);
    s->new_taxname = MemFree (s->new_taxname);
    s = MemFree (s);
  }
  return s;
}


extern ValNodePtr SpecificHostFixListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL)
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = SpecificHostFixFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


extern Boolean ApplyOneSpecificHostFix (SpecificHostFixPtr s)
{
  BioSourcePtr biop = NULL;
  CharPtr      new_spec_host = NULL;
  Boolean      rval = FALSE;

  if (s == NULL || s->feat_or_desc == NULL 
      || StringHasNoText (s->bad_specific_host)
      || StringHasNoText (s->new_taxname)
      || StringHasNoText (s->old_taxname)) return rval;
  biop = GetBioSourceFromValNode (s->feat_or_desc);
  if (biop == NULL) return rval;

  new_spec_host = StringSave (s->bad_specific_host);
  FindReplaceString (&new_spec_host, s->old_taxname, s->new_taxname, TRUE, TRUE);
  if (StringCmp (new_spec_host, s->bad_specific_host) != 0)
  {
    rval = FixSpecificHostForOneBioSource (biop, new_spec_host, s->bad_specific_host);
  }
  new_spec_host = MemFree (new_spec_host);
  return rval;
}

static CharPtr StringIsFirstPartOfItemInList (CharPtr str, ValNodePtr list)
{
  Int4 len;

  if (StringHasNoText (str)) return NULL;
  len = StringLen (str);

  while (list != NULL)
  {
    if (StringNICmp (list->data.ptrvalue, str, len) == 0)
    {
      return list->data.ptrvalue;
    }
    list = list->next;
  }
  return NULL;
}


static ValNodePtr GetFixesForOneSpecificHostValue (SpecificHostCheckPtr p)
{
  CharPtr new_taxname = NULL, old_taxname = NULL;
  OrgRefPtr    request_org, response_org;
  ValNodePtr   biop_vnp, response_vnp, request_vnp;
  SpecificHostFixPtr s;
  ValNodePtr         fix_list = NULL;
  ValNodePtr         failed_requests = NULL;
  CharPtr            whole_name = NULL;

  if (p == NULL) return NULL;

  request_vnp = p->request_list;
  response_vnp = p->response_list;
  
  while (request_vnp != NULL && response_vnp != NULL && old_taxname == NULL)
  {
    request_org = (OrgRefPtr) request_vnp->data.ptrvalue;
    response_org = (OrgRefPtr) response_vnp->data.ptrvalue;
    if (response_org == NULL)
    {
      if (!FindMatchInResponseOrgList (request_org->taxname, p->response_list)) {
        /* no data */
        ValNodeAddPointer (&failed_requests, 0, StringSave (request_org->taxname));
      }
    }
    else if (MatchWithSp (request_org->taxname, response_org->taxname))
    {
      /* found taxname */
      /* is it just a shorter version of an earlier failed request? */
      whole_name = StringIsFirstPartOfItemInList (request_org->taxname, failed_requests);
      if (whole_name != NULL)
      {
        old_taxname = whole_name;
        new_taxname = response_org->taxname; 
      }
      else
      {
        failed_requests = ValNodeFreeData (failed_requests);
        return NULL; 
      }
    }
    else if (!MatchesCommonName (request_org->taxname, response_org->common)
             && !MatchesSynonym (request_org->taxname, response_org))
    {
      old_taxname = request_org->taxname;
      new_taxname = response_org->taxname; 
    }
    request_vnp = request_vnp->next;
    response_vnp = response_vnp->next;
  }

  for (biop_vnp = p->biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {
    s = SpecificHostFixNew (biop_vnp, p->spec_host, old_taxname, new_taxname);
    ValNodeAddPointer (&fix_list, 0, s);
  }
  failed_requests = ValNodeFreeData (failed_requests);
  return fix_list;
}


NLM_EXTERN ValNodePtr Taxon3GetSpecificHostFixesInSeqEntry (SeqEntryPtr sep, Boolean check_single_word)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   req_host_list = NULL, spec_host_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   check_list, check_vnp;
  SpecificHostCheckPtr p;
  ErrSev               level;
  ValNodePtr           fix_list = NULL;
  
  biop_list = GetSpecificHostBioSourceList (sep, FALSE, check_single_word);

  /* get a list of unique specific_host values */
  spec_host_list = GetListOfUniqueSpecificHostValues (biop_list);

  /* now format requests for unique specific_host values */
  FormatSpecificHostRequests (spec_host_list, &request_list, &req_host_list);

  spec_host_list = ValNodeFree (spec_host_list);

  level = ErrSetMessageLevel (SEV_MAX);
  response_list = Taxon3GetOrgRefList (request_list);
  ErrSetMessageLevel (level);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
  }
  else
  {
    /* resort requests so that we can check all responses for the same BioSource together */
    check_list = SortSpecificHostOrgs (req_host_list, request_list, response_list);
    AddBioSourcesToSpecificHostChecklist (biop_list, check_list);  

    /* now look at responses */
    check_vnp = check_list;
    while (check_vnp != NULL)
    {
      p = (SpecificHostCheckPtr) check_vnp->data.ptrvalue;
      ValNodeLink (&fix_list, GetFixesForOneSpecificHostValue (p));
      check_vnp = check_vnp->next;
    }
    check_list = SpecificHostCheckListFree (check_list);
  }

  biop_list = ValNodeFree (biop_list);
  request_list = FreeOrgRefValNodeList (request_list);
  response_list = FreeOrgRefValNodeList (response_list);


  return fix_list;
}


NLM_EXTERN Boolean Taxon3FixSpecificHostInSeqEntry (SeqEntryPtr sep, Boolean check_single_word)
{
  ValNodePtr fix_list, vnp;
  Boolean    rval = FALSE;

  fix_list = Taxon3GetSpecificHostFixesInSeqEntry (sep, check_single_word);
  for (vnp = fix_list; vnp != NULL; vnp = vnp->next)
  {
    rval |= ApplyOneSpecificHostFix (vnp->data.ptrvalue);
  }
  fix_list = SpecificHostFixListFree (fix_list);
  return rval;
}


static void AddBioSourceFeatToList (SeqFeatPtr sfp, Pointer userdata)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;

  ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
}


static void AddBioSourceDescToList (SeqDescrPtr sdp, Pointer userdata)
{

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;

  ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQDESC, sdp);
}


static ValNodePtr GetBioSourceList (SeqEntryPtr sep)
{
  ValNodePtr list = NULL;
  
  VisitFeaturesInSep (sep, &list, AddBioSourceFeatToList);
  VisitDescriptorsInSep (sep, &list, AddBioSourceDescToList);
  return list;
}


static ValNodePtr GetListOfOrganismNames (ValNodePtr biop_list)
{
  ValNodePtr   biop_vnp;
  BioSourcePtr biop;
  ValNodePtr   list = NULL;
  
  /* get a list of unique specific_host values */
  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {
    if (biop_vnp->data.ptrvalue == NULL) continue;
    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL || biop->org == NULL || StringHasNoText (biop->org->taxname)) continue;
    if (!StringAlreadyInValNodeList (biop->org->taxname, list))
    {
      ValNodeAddPointer (&list, 0, biop->org->taxname);
    }
  }
  return list;
}


static void AddBioSourcesToChecklist (ValNodePtr biop_list, ValNodePtr check_list)
{
  ValNodePtr biop_vnp, last_vnp = NULL, stop_search;
  BioSourcePtr biop;
  SpecificHostCheckPtr p;

  if (biop_list == NULL || check_list == NULL) return;

  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {

    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL) continue;

    if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) continue;
    if (last_vnp == NULL)
    {
      last_vnp = check_list;
      stop_search = NULL;
    }
    else
    {
      stop_search = last_vnp;
    }
    p = NULL;
    while (last_vnp != NULL 
           && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
           && StringCmp (p->spec_host, biop->org->taxname) != 0)
    {
      p = NULL;
      last_vnp = last_vnp->next;
    }
    if (p == NULL && stop_search != NULL)
    {
      last_vnp = check_list;
      while (last_vnp != stop_search 
              && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
              && StringCmp (p->spec_host, biop->org->taxname) != 0)
      {
        p = NULL;
        last_vnp = last_vnp->next;
      }
    }

    if (p != NULL)
    {
      ValNodeAddPointer (&(p->biop_list), biop_vnp->choice, biop_vnp->data.ptrvalue);
    }
  }
}


static ValNodePtr GetBioSourcesWithTaxName (CharPtr taxname, ValNodePtr biop_list)
{
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioSourcePtr biop;
  ValNodePtr match_list = NULL, vnp;

  if (StringHasNoText (taxname) || biop_list == NULL) return NULL;

  for (vnp = biop_list; vnp != NULL; vnp = vnp->next) {
    biop = NULL;
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      }
    } else if (vnp->choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      if (sdp != NULL && sdp->choice == Seq_descr_source) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      }
    }
    if (biop != NULL && biop->org != NULL && StringCmp (taxname, biop->org->taxname) == 0) {
      ValNodeAddPointer (&match_list, vnp->choice, vnp->data.ptrvalue);
    }
  }
  return match_list;
}


NLM_EXTERN ValNodePtr GetOrganismTaxLookupFailuresInSeqEntry (SeqEntryPtr sep)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   unique_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   req_vnp, resp_vnp;
  ErrSev               level;
  ValNodePtr           failed_list = NULL, vnp;
  OrgRefPtr            request_org;
  
  biop_list = GetBioSourceList (sep);

  /* get a list of unique specific_host values */
  unique_list = GetListOfOrganismNames (biop_list);

  /* now format requests for unique taxname values */
  for (vnp = unique_list; vnp != NULL; vnp = vnp->next) 
  {
    request_org = OrgRefNew();
    request_org->taxname = StringSave (vnp->data.ptrvalue);
    ValNodeAddPointer (&request_list, 3, request_org);
  }

  unique_list = ValNodeFree (unique_list);

  level = ErrSetMessageLevel (SEV_MAX);
  response_list = Taxon3GetOrgRefList (request_list);
  ErrSetMessageLevel (level);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
  }
  else
  {
    for (req_vnp = request_list, resp_vnp = response_list;
         req_vnp != NULL && resp_vnp != NULL;
         req_vnp = req_vnp->next, resp_vnp = resp_vnp->next)
    {
      if (resp_vnp->data.ptrvalue == NULL)
      {        
        request_org = (OrgRefPtr) req_vnp->data.ptrvalue;
        vnp = GetBioSourcesWithTaxName (request_org->taxname, biop_list);
        if (vnp != NULL) {
          ValNodeAddPointer (&failed_list, 0, StringSave (request_org->taxname));
          ValNodeLink (&failed_list, vnp);
        }
      }
    }
  }

  biop_list = ValNodeFree (biop_list);
  request_list = FreeOrgRefValNodeList (request_list);
  response_list = FreeOrgRefValNodeList (response_list);

  return failed_list;  
}

