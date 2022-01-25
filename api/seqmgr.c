/*  seqmgr.c
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
* File Name:  seqmgr.c
*
* Author:  James Ostell
*   
* Version Creation Date: 9/94
*
* $Revision: 6.54 $
*
* File Description:  Manager for Bioseqs and BioseqSets
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: seqmgr.c,v $
* Revision 6.54  1998/09/30 14:49:37  kans
* set scope for FindAppropriateBioseq, FindFirstLocalBioseq, over and above gather scope, which does not apply to the above calls even from within a gather callback
*
* Revision 6.53  1998/09/29 20:06:07  kans
* FindFirstLocalBioseq and GetOffsetInFirstLocalBioseq to deal with far segments more gracefully than just not indexing the feature
*
* Revision 6.52  1998/09/29 15:07:17  kans
* corrected logic for seqDescFilter, seqFeatFilter, and featDefFilter in explore functions
*
* Revision 6.51  1998/09/23 16:41:07  kans
* added SeqMgrGetDesiredDescriptor
*
* Revision 6.50  1998/09/22 18:17:01  kans
* descriptor index flag now tracked properly, separate from itemID
*
* Revision 6.49  1998/09/22 18:01:25  kans
* had been skipping bssp for lastDescrItemID
*
* Revision 6.48  1998/09/22 16:55:51  kans
* added SeqMgrGetDesiredFeature and position index field
*
* Revision 6.47  1998/09/22 13:11:59  kans
* locationFilter parameter to explore features function
*
* Revision 6.46  1998/09/01 19:25:25  kans
* context parameter in get best protein, get cds/rna given product
*
* Revision 6.45  1998/08/24 18:27:09  kans
* removed solaris -v -fd warnings
*
* Revision 6.44  1998/08/21 21:32:36  kans
* populate ivals array (start/stop pairs) for indexed features
*
* Revision 6.43  1998/08/21 20:18:59  kans
* added SeqMgrExploreSegments, indexing features on segmented bioseq
*
* Revision 6.42  1998/08/19 16:26:48  kans
* MakeReversedSeqIdString called from original location of code, also finished support for biosource feature indexing
*
* Revision 6.41  1998/08/18 21:43:54  kans
* SeqIdWithinBioseq finds appropriate SeqID in bsp->id chain for use with SeqLocAinB, allowing multiple IDs on a protein bioseq
*
* Revision 6.40  1998/08/16 22:36:24  kans
* fixed direct map up from part to segmented bioseq
*
* Revision 6.39  1998/08/14 15:40:37  kans
* SeqMgrMapPartToSegmentedBioseq neede LIBCALL, speeded up function by adding map up on part if fetched
*
* Revision 6.38  1998/08/13 22:31:45  kans
* SeqMgrMapPartToSegmentedBioseq to speed up GetOffsetInBioseq, start of indexing segments, also index biosource by location for binary search (Wheelan)
*
* Revision 6.37  1998/08/12 22:16:29  kans
* sort seg-parts array by SeqId, handle seqloc_int and seqloc_whole
*
* Revision 6.36  1998/08/12 21:25:04  kans
* forgot to free allocated partslist array
*
* Revision 6.35  1998/08/12 21:10:37  kans
* added parts index to speed segmented bioseq mapping
*
* Revision 6.34  1998/07/23 17:30:41  kans
* get overlapping gene was nulling out sfp, then dereferencing
*
* Revision 6.33  1998/07/23 13:08:56  kans
* SeqMgrGetOverlappingGene and SeqMgrGetOverlappingPub take optional context pointer
*
* Revision 6.32  1998/07/23 01:15:19  kans
* minor fix to not return ignored feature
*
* Revision 6.31  1998/07/23 01:11:33  kans
* added SeqMgrGetOverlappingPub, gene overlap works even when gene spans circular origin
*
* Revision 6.30  1998/07/16 22:30:55  kans
* improved gene overlap function
*
* Revision 6.29  1998/07/16 16:56:38  kans
* added parent BioseqSetPtr field to SeqMgrBioseqContext, check most recent bioseq first when indexing features
*
* Revision 6.28  1998/07/06 16:15:17  kans
* fixed typo in explore bioseqs callback
*
* Revision 6.27  1998/07/06 15:57:29  kans
* SeqMgrExploreBioseqs takes entityID or ptr
*
* Revision 6.26  1998/07/06 15:30:15  kans
* scope on index explore, added SeqMgrExploreBioseqs
*
* Revision 6.25  1998/07/02 22:30:44  kans
* process product before indexing by location, ErrPostItem if cannot find bioseq for location
*
* Revision 6.24  1998/07/02 17:52:33  kans
* CreateBioseqExtraBlock was not being called for protein bioseq to link back to CDS, which was seen first
*
* Revision 6.23  1998/07/01 19:13:16  kans
* SMFeatBlock.data is allocated array of reasonable size
*
* Revision 6.22  1998/06/30 22:08:02  kans
* SeqMgrFeaturesAreIndexed takes entityID, returns time_t time stamp of latest indexing
*
* Revision 6.21  1998/06/30 14:28:00  kans
* changed GetSeqFeat, which collided with asn2ff4 for some linkers
*
* Revision 6.20  1998/06/30 14:20:16  kans
* changes to heap sort order to put genes first, then rnas, if ranges are equal
*
* Revision 6.19  1998/06/30 12:56:45  kans
* code fixes, public functions moved to explore.h
*
* Revision 6.18  1998/06/29 23:37:37  kans
* added context structure for all explores, index every bioseq in an entity
*
* Revision 6.17  1998/06/29 03:06:41  kans
* fixed two conditionals in ProcessFeatureProducts
*
* Revision 6.16  1998/06/29 02:29:44  kans
* new context get descriptor and feature functions now working, get gene not always working
*
* Revision 6.15  1998/06/29 01:33:27  kans
* added SeqMgrGetNextDescriptor and SeqMgrGetNextFeature
*
* Revision 6.14  1998/06/29 00:24:00  kans
* several changes to new indexing functions
*
* Revision 6.13  1998/06/28 03:44:18  kans
* omdp->parentptr is bssp, not omdp, so use ObjMgrFindByData to get higher descriptors
*
* Revision 6.12  1998/06/28 03:15:09  kans
* missing break statement caused features to be ignored
*
* Revision 6.11  1998/06/28 02:38:15  kans
* simplified filters, finished best gene, explore functions
*
* Revision 6.10  1998/06/27 22:23:49  kans
* improvements and further implementation of new indexing, exploration functions
*
* Revision 6.9  1998/06/27 00:03:45  kans
* fix feature heap sort, post increment feature insertion index, look for best prot
*   when setting cds back pointer, and merge descriptor count and feature collect callbacks
*
* Revision 6.8  1998/06/26 22:36:24  kans
* initial work on tracking sorted features, and cds and prot links, for rapid collection
*
* Revision 6.7  1998/05/01 16:13:13  kans
* caching of gi with NULL seqID allowed with protection against calling SeqIdDup
*
* Revision 6.6  1998/04/20 22:38:08  kans
* should prevent caching of gi with NULL seqID
*
* Revision 6.5  1998/04/08 16:52:08  kans
* casts to ValNodeLen calls
*
* Revision 6.4  1998/03/30 21:02:25  ostell
* removed check for parenttype != 0 on call to ObjMgrConnect in SeqMgrLinkSeqEntry
*   so that disconnects from sets would work as well as connects
*
* Revision 6.3  1997/11/19 22:14:42  ostell
* added support for multithreaded programs
*
* Revision 6.2  1997/09/25 18:20:14  tatiana
* fixing -1 bug for gaps in CountGapsInDeltaSeq
*
* Revision 6.1  1997/09/11 15:55:40  ostell
* Added support for SetColor messages
*
* Revision 6.0  1997/08/25 18:07:06  madden
* Revision changed to 6.0
*
* Revision 5.20  1997/07/31 16:06:49  kans
* BioseqLockById clears scope if first call to BioseqFindFunc fails, tries again
*
* Revision 5.19  1997/07/30 19:44:46  kans
* bug fix by Serge Bazhin
*
* Revision 5.18  1997/07/28 13:29:41  ostell
* Moved GetUniGeneIDForSeqId() to seqmgr.c
*
* Revision 5.17  1997/07/15 17:37:43  ostell
* fixed problems with duplicate Bioseqs in BioseqFindFunc
*
* Revision 5.16  1997/07/09 21:11:53  ostell
* added support for indexed seqid lookups of bioseqs
*
* Revision 5.15  1997/06/19 18:38:43  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.14  1997/06/17 17:59:09  kans
* GetSeqIdForGI now uses binary search in cache
*
* Revision 5.13  1997/06/17 16:33:57  kans
* first pass at cache in SeqIdForGi
*
* Revision 5.12  1997/03/26 14:01:37  ostell
* removed OMUserData from new copy of Bioseq when uncaching to stop
* memory leak
*
 * Revision 5.11  1997/02/24  21:46:17  ostell
 * in BioseqFindFunc, when checking with scope a failure occurs, now it
 *   checks again without scope.
 *
 * Revision 5.10  1997/01/23  22:38:21  ostell
 * added missing newline at end of file (sigh)
 *
 * Revision 5.9  1997/01/23  22:37:14  ostell
 * minor change to seqmgr.h for new indexing
 *
 * Revision 5.7  1997/01/08  22:48:50  tatiana
 * buf and buflen arguments added to CountGapsInDeltaSeq()
 *
 * Revision 5.6  1996/08/22  14:50:05  ostell
 * initialized static arrays in BioseqFindFunc
 *
 * Revision 5.5  1996/08/21  13:33:33  ostell
 * added cachig to BioseqFindFunc
 *
 * Revision 5.4  1996/08/06  19:56:03  kans
 * for SEQLOC_WHOLE, must call SeqIdFindBest on bsp->id
 *
 * Revision 5.3  1996/08/05  15:57:26  chappey
 * in BioseqReloadFunc, the OMUserDataPtr is passed to the new
 * ObjMgrDataPtr, and is not deleted anymore.
 *
 * Revision 5.2  1996/07/25  02:32:26  ostell
 * added CountGapsInDeltaSeq()
 *
 * Revision 5.1  1996/07/19  22:13:13  ostell
 * added SpreadGapsInDeltaSeq()
 *
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 4.7  1996/03/19  19:05:17  kans
 * SeqEntrySetScope now returns old scope, not new scope
 *
 * Revision 4.6  1996/01/23  14:44:38  kans
 * added Pointer casts to MemSet
 *
 * Revision 4.5  1995/12/22  14:43:59  ostell
 * added reload code to BioseqLockById
 * break out relad from cache code to be used as part of gather locking
 * with BioseqReload
 *
 * Revision 4.4  1995/12/09  23:12:41  kans
 * SeqEntryFind now can deal with a Seq-Submit ultimate parent
 *
 * Revision 4.3  1995/12/04  21:40:05  ostell
 * added GetSeqIdForGI() and GetGIForSeqId()
 *
 * Revision 4.2  1995/10/03  15:50:37  ostell
 * added support for selection by region.. now fully implemented
 *
 * Revision 4.1  1995/09/30  03:38:31  ostell
 * Changed ObjMgrMessage functions to pass a structure
 * Added support for selecting regions
 * Added ability to remove entity when no more views on it
 *
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 1.16  1995/05/15  21:46:05  ostell
 * added Log line
 *
*
*
*
* ==========================================================================
*/

/** for ErrPostEx() ****/

static char *this_module = "ncbiapi";
#define THIS_MODULE this_module
static char *this_file = __FILE__;
#define THIS_FILE this_file

/**********************/

#include <explore.h>       /* new public functions prototyped here */
#include <seqmgr.h>		   /* the interface */
#include <sequtil.h>       /* CLEAN THIS UP LATER? */
#include <gather.h>
#include <subutil.h>
#include <ncbithr.h>
#include <objfdef.h>
#include <sqnutils.h>

/*****************************************************************************
*
*   Bioseq Management
*
*****************************************************************************/

static BioseqPtr LIBCALLBACK BSFetchFunc PROTO((SeqIdPtr sid, Uint1 ld_type));
static BioseqPtr NEAR BioseqFindFunc PROTO((SeqIdPtr sid, Boolean reload_from_cache));
static Boolean NEAR SeqMgrGenericSelect PROTO((SeqLocPtr region, Int2 type,
                                             Uint1Ptr rgb));
static BioseqPtr NEAR BioseqReloadFunc PROTO((SeqIdPtr sid, ObjMgrDataPtr omdp));

static Boolean NEAR SeqMgrProcessNonIndexedBioseq PROTO((void));
static Boolean NEAR SeqMgrAddIndexElement PROTO((SeqMgrPtr smp, BioseqPtr bsp, CharPtr buf));
static void NEAR RevStringUpper PROTO((CharPtr str));
static BSFetchTop NEAR SeqMgrGetFetchTop (void);


/*****************************************************************************
*
*   Return the current SeqMgr
*   	SeqMgrGet is obsolete
*       SeqMgrReadLock, ReadUnlock, WriteLock, WriteUnlock are thread safe
*
*****************************************************************************/
static TNlmMutex smp_mutex = NULL;
static SeqMgrPtr global_smp = NULL;
static TNlmRWlock smp_RWlock = NULL;
static TNlmRWlock sgi_RWlock = NULL;

/*****************************************************************************
*
*   Return the current SeqMgr
*   	Initialize if not done already
*       This function will become obsolete
*
*****************************************************************************/
NLM_EXTERN SeqMgrPtr LIBCALL SeqMgrGet (void)
{
	Int4 ret;
	SeqMgrPtr smp;

	if (global_smp != NULL)
		return global_smp;

	ret = NlmMutexLockEx(&smp_mutex);  /* protect this section */
	if (ret)  /* error */
	{
		ErrPostEx(SEV_FATAL,0,0,"SeqMgrGet failed [%ld]", (long)ret);
		return NULL;
	}

	if (global_smp == NULL)  /* check again after mutex */
	{
	                             /*** have to initialize it **/
		smp = (SeqMgrPtr) MemNew (sizeof(SeqMgr));
		smp->bsfetch = BSFetchFunc;  /* BioseqFetch default */
		smp->fetch_on_lock = TRUE;	 /* fetch when locking */
		smp_RWlock = NlmRWinit();  /* initialize RW lock */
		sgi_RWlock = NlmRWinit();  /* initialize RW lock */
		global_smp = smp;       /* do this last for mutex safety */
	}

	NlmMutexUnlock(smp_mutex);

	return global_smp;
}

/*****************************************************************************
*
*   SeqMgrReadLock()
*   	Initialize if not done already
*       A thread can have only one read or write lock at a time
*       Many threads can have read locks
*       Only one thread can have a write lock
*       No other threads may have read locks if a write lock is granted
*       If another thread holds a write lock, this call blocks until write
*          is unlocked.
*
*****************************************************************************/
NLM_EXTERN SeqMgrPtr LIBCALL SeqMgrReadLock (void)
{
	SeqMgrPtr smp;
	Int4 ret;

	smp = SeqMgrGet();  /* ensure initialization */

	ret = NlmRWrdlock(smp_RWlock);
	if (ret != 0)
	{
		ErrPostEx(SEV_ERROR,0,0,"SeqMgrReadLock: RWrdlock error [%ld]",
			(long)ret);
		return NULL;
	}
	return smp;
}

/*****************************************************************************
*
*   SeqMgrWriteLock
*   	Initialize if not done already
*       A thread can have only one read or write lock at a time
*       Many threads can have read locks
*       Only one thread can have a write lock
*       No other threads may have read locks if a write lock is granted
*       If another thread holds a read or write lock, this call blocks until write
*          is unlocked.
*
*****************************************************************************/
NLM_EXTERN SeqMgrPtr LIBCALL SeqMgrWriteLock (void)
{
	SeqMgrPtr smp;
	Int4 ret;

	smp = SeqMgrGet();  /* ensure initialization */

	ret = NlmRWwrlock(smp_RWlock);
	if (ret != 0)
	{
		ErrPostEx(SEV_ERROR,0,0,"SeqMgrWriteLock: RWwrlock error [%ld]",
			(long)ret);
		return NULL;
	}
	smp->is_write_locked = TRUE;
	return smp;
}


/*****************************************************************************
*
*  SeqMgrUnlock()
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrUnlock (void)
{
	SeqMgrPtr smp;
	Int4 ret;

	smp = SeqMgrGet();  /* ensure initialization */

	ret = NlmRWunlock(smp_RWlock);
	if (ret != 0)
	{
		ErrPostEx(SEV_ERROR,0,0,"SeqMgrUnlock: RWunlock error [%ld]",
			(long)ret);
		return FALSE;
	}
	smp->is_write_locked = FALSE;  /* can't be write locked */
	return TRUE;
}

/****************************************************************************
*
*  RevStringUpper(str)
*    Up cases and reverses string
*      to get different parts early for SeqId StringCmps
*
*****************************************************************************/
static void NEAR RevStringUpper (CharPtr str)
{
	CharPtr nd;
	Char tmp;

		if (str == NULL)
			return;
    nd = str;
	while (*nd != '\0')
		nd++;
	nd--;

	while (nd > str)
	{
		tmp = TO_UPPER(*nd);
		*nd = TO_UPPER(*str);
		*str = tmp;
		nd--; str++;
	}

	if (nd == str)
		*nd = TO_UPPER(*nd);
	return;
}

static Boolean MakeReversedSeqIdString (SeqIdPtr sid, CharPtr buf, size_t len)

{
  Uint1         oldchoice;
  CharPtr       tmp;
  TextSeqIdPtr  tsip;

  if (sid == NULL || buf == NULL || len < 1) return FALSE;
  oldchoice = 0;
  switch (sid->choice) {
    case SEQID_GI:
      sprintf (buf, "%ld", (long)(sid->data.ptrvalue));
      break;
    case SEQID_EMBL:
    case SEQID_DDBJ:
      oldchoice = sid->choice;
      sid->choice = SEQID_GENBANK;
    case SEQID_GENBANK:
    case SEQID_PIR:
    case SEQID_OTHER:
    case SEQID_SWISSPROT:
    case SEQID_PRF:
      tsip = (TextSeqIdPtr) (sid->data.ptrvalue);
      if (tsip->accession != NULL) {
        tmp = tsip->name;
        tsip->name = NULL;
        SeqIdWrite (sid, buf, PRINTID_FASTA_SHORT, len);
        tsip->name = tmp;
      } else {
        SeqIdWrite (sid, buf, PRINTID_FASTA_SHORT, len);
      }
      if (oldchoice)
        sid->choice = oldchoice;
      break;
    default:
      SeqIdWrite (sid, buf, PRINTID_FASTA_SHORT, len);
      break;
  }
  RevStringUpper (buf);
  return TRUE;
}

/*****************************************************************************
*
*   SeqEntrySetScope(sep)
*   	scopes global seqentry searches to sep
*       setting sep=NULL, opens scope to all seqentries in memory
*       returns the current scope
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr LIBCALL SeqEntrySetScope(SeqEntryPtr sep)
{
	SeqEntryPtr curr = NULL;
	SeqMgrPtr smp;
	Int2 i, j;
	SMScopePtr smsp;
	TNlmThread thr;
	Boolean found;

	smp = SeqMgrWriteLock();
	if (smp == NULL) goto ret;
	thr = NlmThreadSelf();
	found = FALSE;
	for (i = 0, smsp = smp->scope; i < smp->num_scope; i++, smsp++)
	{
		if (NlmThreadCompare(thr, smsp->thr))
		{
			curr = smsp->SEscope;
			smsp->SEscope = sep;
			if (sep == NULL)  /* removing one? */
			{
				smp->num_scope--;
				j = smp->num_scope - i;  /* number to move */
				if (j)  /* not last one */
					MemCopy(smsp, (smsp+1), (size_t)(j * sizeof(SMScope)));
			}
			goto ret;    /* all done */
		}
	}

	              /* thread not on list */
	if (sep == NULL)
		goto ret;       /* nothing to do */

	i = smp->num_scope;
	j = smp->total_scope;
	if (j == i)  /* need more room */
	{
		j += 20;   /* new size */
		smsp = smp->scope;
		smp->scope = MemNew((size_t)(j * sizeof(SMScope)));
		MemCopy(smp->scope, smsp, (size_t)(i * sizeof(SMScope)));
		smp->total_scope = j;
		MemFree(smsp);
	}

	smp->scope[i].thr = thr;
	smp->scope[i].SEscope = sep;
	smp->num_scope++;

ret: SeqMgrUnlock();
	return curr;
}

/*****************************************************************************
*
*   SeqEntryGetScope(sep)
*       returns the current scope or NULL if none set
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr LIBCALL SeqEntryGetScope(void)
{
	SeqMgrPtr smp;
	SeqEntryPtr scope = NULL;
	Int2 i;
	SMScopePtr smsp;
	TNlmThread thr;

	smp = SeqMgrReadLock();
	if (smp == NULL) return FALSE;
	thr = NlmThreadSelf();
	for (i = 0, smsp = smp->scope; i < smp->num_scope; i++, smsp++)
	{
		if (NlmThreadCompare(thr, smsp->thr))
		{
			scope = smsp->SEscope;
			break;
		}
	}
	SeqMgrUnlock();
	return scope;
}

/*****************************************************************************
*
*   BioseqFind(SeqIdPtr)
*   	Just checks in object loaded memory
*       Will also restore a Bioseq that has been cached out
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqFind (SeqIdPtr sid)
{
	return BioseqFindFunc(sid, TRUE);
}

/*****************************************************************************
*
*   BioseqFindCore(sid)
*   	Finds a Bioseq in memory based on SeqId when only "core" elements needed
*   	Will NOT restore a Bioseq that has been cached out by SeqMgr
*       This function is for use ONLY by functions that only need the parts
*         of the Bioseq left when cached out. This includes the SeqId chain,
*         and non-pointer components of the Bioseq.
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqFindCore (SeqIdPtr sip)
{
	return BioseqFindFunc(sip, FALSE);
}

/*****************************************************************************
*
*   BioseqFindEntity(sid, itemIDptr)
*   	Finds a Bioseq in memory based on SeqId
*   	Will NOT restore a Bioseq that has been cached out by SeqMgr
*       returns EntityID if found, otherwise 0
*       itemIDptr is set to the value for itemID in ObjMgr functions
*       itemtype is OBJ_BIOSEQ of course
*
*****************************************************************************/
NLM_EXTERN Uint2 LIBCALL BioseqFindEntity (SeqIdPtr sip, Uint2Ptr itemIDptr)
{
	BioseqPtr bsp;
	Uint2 entityID = 0;

	*itemIDptr = 0;
	bsp = BioseqFindCore(sip);
	if (bsp == NULL) return entityID;  /* not found */
	entityID = ObjMgrGetEntityIDForPointer((Pointer)bsp);
	if (! entityID)
		return entityID;

	*itemIDptr = GatherItemIDByData(entityID, OBJ_BIOSEQ, (Pointer)bsp);
	return entityID;
}

/********************************************************************************
*
*   BioseqReload (omdp, lockit)
*     reloads the cached SeqEntry at top of omdp
*     if (lockit) locks the record
*
*********************************************************************************/

NLM_EXTERN ObjMgrDataPtr LIBCALL BioseqReload(ObjMgrDataPtr omdp, Boolean lockit)
{
	BioseqPtr bsp = NULL;
	ObjMgrDataPtr retval = NULL;
	Int2 j;
	ObjMgrPtr omp;

	if (omdp == NULL) return retval;
	if (! ((omdp->datatype == OBJ_BIOSEQ) || (omdp->datatype == OBJ_BIOSEQSET)))
		return retval;
	if (omdp->parentptr != NULL)
	{
		omp = ObjMgrReadLock();
		omdp = ObjMgrFindTop(omp, omdp);
		ObjMgrUnlock();
		if (omdp == NULL)
			return retval;
	}

	if (omdp->tempload == TL_CACHED)   /* only need to reload if cached */
	{
		bsp = BioseqReloadFunc (NULL, omdp);
		if (bsp == NULL)
			return retval;
		omp = ObjMgrReadLock();
		j = ObjMgrLookup(omp, (Pointer)bsp);
		omdp = ObjMgrFindTop(omp, omp->datalist[j]);
		ObjMgrUnlock();
	}
	 
	if (lockit)
	{
		ObjMgrLock(omdp->datatype, omdp->dataptr, TRUE);
	}

	return omdp;
}

static BSFetchTop NEAR SeqMgrGetFetchTop (void)
{
	SeqMgrPtr smp;
	BSFetchTop bsftp=NULL;

	smp = SeqMgrReadLock();
	if (smp == NULL) return bsftp;
	bsftp = smp->bsfetch;
	SeqMgrUnlock();
	return bsftp;
}
	
static BioseqPtr NEAR BioseqReloadFunc (SeqIdPtr sid, ObjMgrDataPtr omdp)
{
    Int2 j;
	ObjMgrDataPtr oldomdp;
	OMUserDataPtr omudp, next;
	ObjMgrProcPtr ompp;
	OMProcControl ompc;
	BioseqPtr bsp= NULL;
	Int2 ret;
	ObjMgrPtr omp;
	BSFetchTop bsftp=NULL;

	ompp = NULL;
	omp = ObjMgrReadLock();
	for (omudp = omdp->userdata; omudp != NULL; omudp = omudp->next)
	{
		if (omudp->proctype == OMPROC_FETCH)  /* caching function */
		{
			ompp = ObjMgrProcFind(omp, omudp->procid, NULL, 0);
			if (ompp != NULL)
				break;
		}
	}
	ObjMgrUnlock();

	if (ompp == NULL)
		return bsp;
	if (ompp->outputtype != OBJ_BIOSEQ)
		return bsp;

	oldomdp = omdp;
	omdp = NULL;
	bsftp = SeqMgrGetFetchTop();
	if (bsftp != NULL)
	{
		if (ompp != NULL)	/* fetch proc left a signal */
		{                                 /* rerun fetch */
			MemSet((Pointer)(&ompc), 0, sizeof(OMProcControl));
			ompc.input_data = sid;
			ompc.input_entityID = oldomdp->EntityID;
			ompc.proc = ompp;
			ret = (* (ompp->func))((Pointer)&ompc);
			switch (ret)
			{
				case OM_MSG_RET_ERROR:
					ErrShow();
					break;
				case OM_MSG_RET_DEL:
					break;
				case OM_MSG_RET_OK:
					break;
				case OM_MSG_RET_DONE:
					omp = ObjMgrWriteLock();
					ObjMgrSetTempLoad (omp, ompc.output_data);
					ObjMgrUnlock();
					bsp = (BioseqPtr)(ompc.output_data);
					break;
				default:
					break;
			}
		}
		
		if (bsp == NULL)  /* nope, try regular fetch */
		{
			bsp = (*(bsftp))(sid, BSFETCH_TEMP);
		}

		if (bsp != NULL)
		{
			omp = ObjMgrReadLock();
			j = ObjMgrLookup(omp, (Pointer)bsp);
			omdp = ObjMgrFindTop(omp, omp->datalist[j]);
			ObjMgrUnlock();
			omdp->EntityID = oldomdp->EntityID;
			oldomdp->EntityID = 0;

			omudp = omdp->userdata;
			while (omudp != NULL)
			{
				next = omudp->next;
				if (omudp->freefunc != NULL)
                                 (*(omudp->freefunc))(omudp->userdata.ptrvalue);
				MemFree(omudp);
				omudp = next;
			}
            omdp->userdata = oldomdp->userdata;
            oldomdp->userdata = NULL;

			if (oldomdp->choice != NULL)
				SeqEntryFree(oldomdp->choice);
			else
			{
				switch(oldomdp->datatype)
				{
					case OBJ_BIOSEQ:
						BioseqFree((BioseqPtr)(oldomdp->dataptr));
						break;
					case OBJ_BIOSEQSET:
						BioseqSetFree((BioseqSetPtr)(oldomdp->dataptr));
						break;
					default:
						ErrPostEx(SEV_ERROR,0,0,"BioseqFindFunc: delete unknown type [%d]",
							(int)(oldomdp->datatype));
						break;
				}
			}
		}
	}
	return bsp;
}
/** static func used internally **/

/*******************************************
*
*  WARNING: if you change BIOSEQ_CACHE_NUM, you have to change the
*   number of NULL in the initialization of the 2 static pointer arrays
*   below
*
*******************************************/
/* nb: this cache is cleared in SeqMgrDeleteFromBioseqIndex() */
#define BIOSEQ_CACHE_NUM 3
static SeqEntryPtr se_cache[BIOSEQ_CACHE_NUM] = {
	NULL, NULL, NULL};   /* for a few platforms */
static ObjMgrDataPtr omdp_cache[BIOSEQ_CACHE_NUM] = {
	NULL, NULL, NULL};   /* for a few platforms */
static TNlmMutex smp_cache_mutex = NULL;

static BioseqPtr NEAR BioseqFindFunc (SeqIdPtr sid, Boolean reload_from_cache)
{
    Int4 i, j, num, imin, imax, ret;
	SeqIdIndexElementPtr PNTR sipp;
	CharPtr tmp;
	Char buf[80];
	Boolean do_return;
	SeqMgrPtr smp;
	ObjMgrPtr omp;
	ObjMgrDataPtr omdp;
	BioseqPtr bsp = NULL, tbsp;
	SeqEntryPtr scope;

	if (sid == NULL)
		return NULL;

	ret = NlmMutexLockEx(&smp_cache_mutex);  /* protect this section */
	if (ret)  /* error */
	{
		ErrPostEx(SEV_FATAL,0,0,"BioseqFindFunc cache mutex failed [%ld]", (long)ret);
		return NULL;
	}

	do_return = FALSE;
	scope = SeqEntryGetScope();	   /* first check the cache */
	for (i = 0; i < BIOSEQ_CACHE_NUM; i++)
	{
		if (omdp_cache[i] == NULL)
			break;
		omdp = omdp_cache[i];
		if (omdp->datatype == OBJ_BIOSEQ)
		{
			if ((scope == NULL) || (scope == se_cache[i]))
			{
				bsp = (BioseqPtr)(omdp->dataptr);

				if (BioseqMatch(bsp, sid))
				{
					for (j = i; j > 0; j--)  /* shift to top of cache */
					{
						omdp_cache[j] = omdp_cache[j-1];
						se_cache[j] = se_cache[j-1];
					}
					omdp_cache[0] = omdp;
					se_cache[0] = scope;

					if (! reload_from_cache)
					{
						do_return = TRUE;
						goto done_cache;
					}

					omp = ObjMgrReadLock();
					omdp = ObjMgrFindTop(omp, omdp);
					ObjMgrUnlock();
					if (omdp->tempload != TL_CACHED)
					{
						do_return = TRUE;
						goto done_cache;
					}

					bsp = BioseqReloadFunc(sid, omdp);

					if (bsp == NULL)
					{
						
                        ErrPostEx(SEV_ERROR,0,0,"BioseqFindFunc: couldn't uncache");
					}
					do_return = TRUE;
					goto done_cache;
				}
			}
		}
	}
done_cache:
	NlmMutexUnlock(smp_cache_mutex);
	if (do_return)  /* all done */
	{
		return bsp;
	}

	bsp = NULL; /* resetting it */

	SeqMgrProcessNonIndexedBioseq();	/* make sure all are indexed */

		/* stringify as in SeqMgrAdd */

	MakeReversedSeqIdString (sid, buf, 79); /* common function to make id, call RevStringUpper */

	/*
	oldchoice = 0;
	switch (sid->choice)
	{
	case SEQID_GI:
		sprintf(buf, "%ld", (long)(sid->data.ptrvalue));
		break;
	case SEQID_EMBL:
	case SEQID_DDBJ:
		oldchoice = sid->choice;
		sid->choice = SEQID_GENBANK;
	case SEQID_GENBANK:
	case SEQID_PIR:
	case SEQID_OTHER:
	case SEQID_SWISSPROT:
	case SEQID_PRF:
		tsip = (TextSeqIdPtr)(sid->data.ptrvalue);
		if (tsip->accession != NULL)
		{
			tmp = tsip->name;
			tsip->name = NULL;
			SeqIdWrite(sid, buf, PRINTID_FASTA_SHORT, 79);
			tsip->name = tmp;
		}
		else
 			SeqIdWrite(sid, buf, PRINTID_FASTA_SHORT, 79);
		if (oldchoice)
			sid->choice = oldchoice;
		break;
	default:
  		SeqIdWrite(sid, buf, PRINTID_FASTA_SHORT, 79);
		break;
	}

	RevStringUpper(buf);
	*/

	imin = 0;
	smp = SeqMgrReadLock();
	imax = smp->BioseqIndexCnt - 1;
	sipp = smp->BioseqIndex;

	num = -1;

	while (imax >= imin)
	{
		i = (imax + imin)/2;
		tmp = sipp[i]->str;
		if ((j = StringCmp(tmp, buf)) > 0)
			imax = i - 1;
		else if (j < 0)
			imin = i + 1;
		else
		{
			num = i;
			break;
		}
	}

	if (num < 0)  /* couldn't find it */
	{
		/*
		Message(MSG_ERROR, "[1] Couldn't find [%s]", buf);
		*/
		SeqMgrUnlock();
		return NULL;
	}


	if (scope != NULL)	/* check in scope */
	{
		tbsp = (BioseqPtr)(sipp[num]->omdp->dataptr);
		if (ObjMgrIsChild(scope->data.ptrvalue, tbsp))
		{
			bsp = tbsp;
			omdp = sipp[num]->omdp;
		}
		else
		{                  /* not in scope, could be duplicate SeqId */
			i = num-1;
			while ((i >= 0) && (bsp == NULL) && (! StringCmp(sipp[i]->str, buf)))  /* back up */
			{
			   tbsp = (BioseqPtr)(sipp[i]->omdp->dataptr);
			   if (ObjMgrIsChild(scope->data.ptrvalue, tbsp))
			   {
				   bsp = tbsp;
					omdp = sipp[i]->omdp;
			   }
			   i--;
			}
			i = num + 1;
			imax = smp->BioseqIndexCnt - 1;
			while ((bsp == NULL) && (i <= imax) && (! StringCmp(sipp[i]->str, buf)))
			{
			   tbsp = (BioseqPtr)(sipp[i]->omdp->dataptr);
			   if (ObjMgrIsChild(scope->data.ptrvalue, tbsp))
			   {
				   bsp = tbsp;
					omdp = sipp[i]->omdp;
			   }
			   i++;
			}
		}
	}
	else  /* no scope set */
	{
		omdp = sipp[num]->omdp;
		bsp = (BioseqPtr)(omdp->dataptr);
	}

	SeqMgrUnlock();

	if (bsp == NULL)   /* not found */
	{
		/*
		Message(MSG_ERROR, "[2] Couldn't find [%s]", buf);
		*/
		return bsp;
	}

	ret = NlmMutexLockEx(&smp_cache_mutex);  /* protect this section */
	if (ret)  /* error */
	{
		ErrPostEx(SEV_FATAL,0,0,"BioseqFindFunc2 cache mutex failed [%ld]", (long)ret);
		return NULL;
	}

	for (j = (BIOSEQ_CACHE_NUM - 1); j > 0; j--)  /* shift to top of cache */
	{
		omdp_cache[j] = omdp_cache[j-1];
		se_cache[j] = se_cache[j-1];
	}
	omdp_cache[0] = omdp;
	se_cache[0] = scope;

	NlmMutexUnlock(smp_cache_mutex);

	if (! reload_from_cache)
		return bsp;

	omp = ObjMgrReadLock();
	omdp = ObjMgrFindTop(omp, omdp);
	ObjMgrUnlock();
	if (omdp == NULL)
	{
		bsp = NULL;
		goto ret;
	}
	if (omdp->tempload == TL_CACHED)
		bsp = BioseqReloadFunc(sid, omdp);
ret:
	return bsp;
}

/*****************************************************************************
*
*   SeqMgrFreeCache()
*   	frees all cached SeqEntrys
*   	returns FALSE if any errors occurred
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrFreeCache(void)
{
	return ObjMgrFreeCache(OBJ_SEQENTRY);
}

/*****************************************************************************
*
*   BioseqLockById(SeqIdPtr)
*   	Finds the Bioseq and locks it
*       Makes sure appropriate BioseqContent is present
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqLockById (SeqIdPtr sid)
{
	BioseqPtr bsp = NULL;
	SeqMgrPtr smp;
	SeqEntryPtr oldscope = NULL;
	BSFetchTop bsftp;
	Boolean fetch_on_lock;

	if (sid == NULL) return bsp;

	bsp = BioseqFindFunc(sid, TRUE);
	if (bsp == NULL)
	{
		smp = SeqMgrReadLock();
		if (smp == NULL) return bsp;
		fetch_on_lock = smp->fetch_on_lock;
		bsftp = smp->bsfetch;
		SeqMgrUnlock();

		if (fetch_on_lock)
		{
			oldscope = SeqEntrySetScope (NULL);
			if (oldscope != NULL) {
				bsp = BioseqFindFunc(sid, TRUE);
				SeqEntrySetScope (oldscope);
			}
			if (bsp == NULL && bsftp != NULL)
	            bsp = (*(bsftp))(sid, BSFETCH_TEMP);
		}
	}

	if (bsp == NULL) return NULL;

	ObjMgrLock(OBJ_BIOSEQ, (Pointer)bsp, TRUE);
	return bsp;
}

/*****************************************************************************
*
*   BioseqUnlockById(SeqIdPtr sip)
*   	Frees a Bioseq to be dumped from memory if necessary
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqUnlockById (SeqIdPtr sip)
{
	BioseqPtr bsp;

	if (sip == NULL) return FALSE;

	bsp = BioseqFindFunc(sip, FALSE);
	if (bsp == NULL)
		return FALSE;

	ObjMgrLock(OBJ_BIOSEQ, (Pointer)bsp, FALSE);
	return TRUE;
}

/*****************************************************************************
*
*   BioseqLock(BioseqPtr)
*   	Locks a Bioseq
*       Any cached data is returned to memory
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqLock (BioseqPtr bsp)
{
	if (bsp == NULL) return NULL;

	ObjMgrLock(OBJ_BIOSEQ, (Pointer)bsp, TRUE);

	return bsp;
}

/*****************************************************************************
*
*   BioseqUnlock(BioseqPtr)
*   	Frees a Bioseq to be dumped from memory if necessary
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqUnlock (BioseqPtr bsp)
{
	if (bsp == NULL) return FALSE;

	if (ObjMgrLock(OBJ_BIOSEQ, (Pointer)bsp, FALSE) >= 0)
		return TRUE;
	else
		return FALSE;
}

/*****************************************************************************
*
*   BioseqFetch(SeqIdPtr, flag)
*   	loads bioseq into memory if possible
*   	first trys LocalLoad
*       they trys EntrezLoad
*
*****************************************************************************/
static BioseqPtr LIBCALLBACK BSFetchFunc (SeqIdPtr sid, Uint1 ld_type)
{
	BioseqPtr bsp = NULL;
	ObjMgrProcPtr ompp;
	OMProcControl ompc;
	Int2 ret;
	ObjMgrPtr omp;

	ompp = NULL;
	while ((ompp = ObjMgrProcFindNext(NULL, OMPROC_FETCH, OBJ_SEQID, OBJ_BIOSEQ, ompp)) != NULL)
	{
		MemSet((Pointer)(&ompc), 0, sizeof(OMProcControl));
		ompc.input_data = sid;
		ompc.proc = ompp;
		ret = (* (ompp->func))((Pointer)&ompc);
		switch (ret)
		{
			case OM_MSG_RET_ERROR:
				ErrShow();
				break;
			case OM_MSG_RET_DEL:
				break;
			case OM_MSG_RET_OK:
				break;
			case OM_MSG_RET_DONE:
				if (ld_type == BSFETCH_TEMP)
				{
					omp = ObjMgrWriteLock();
					ObjMgrSetTempLoad (omp, ompc.output_data);
					ObjMgrUnlock();
				}
				bsp = (BioseqPtr)(ompc.output_data);
				break;
			default:
				break;
		}
		if (bsp != NULL)  /* got one */
			break;
	}

	return bsp;
}


NLM_EXTERN BioseqPtr LIBCALL BioseqFetch (SeqIdPtr sid, Uint1 ld_type)
{
	BSFetchTop bsftp;
	BioseqPtr bsp;

	bsp = BioseqFindFunc(sid, TRUE);
	if (bsp != NULL) return bsp;
	
	bsftp = SeqMgrGetFetchTop();
	if (bsftp == NULL) return NULL;

	return (*(bsftp))(sid, ld_type);
}

/*****************************************************************************
*
*   GetGIForSeqId(SeqIdPtr)
*     returns the GI for a SeqId
*     returns 0 if can't find it
*
*****************************************************************************/
NLM_EXTERN Int4 LIBCALL GetGIForSeqId (SeqIdPtr sid)
{
	BioseqPtr bsp = NULL;
	ObjMgrProcPtr ompp;
	OMProcControl ompc;
	Int2 ret;
	SeqIdPtr sip;
	Int4 gi=0;


	if (sid == NULL)
		return gi;

	if (sid->choice == SEQID_GI)
		return sid->data.intvalue;

	bsp = BioseqFindCore(sid);
	if (bsp != NULL)
	{
		for (sip = bsp->id; sip != NULL; sip = sip->next)
		{
			if (sip->choice == SEQID_GI)
				return sip->data.intvalue;
		}
	}


	ompp = NULL;
	while ((ompp = ObjMgrProcFindNext(NULL, OMPROC_FETCH, OBJ_SEQID, OBJ_SEQID, ompp)) != NULL)
	{
		if ((ompp->subinputtype == 0) && (ompp->suboutputtype == SEQID_GI))
		{
			MemSet((Pointer)(&ompc), 0, sizeof(OMProcControl));
			ompc.input_data = sid;
			ompc.proc = ompp;
			ret = (* (ompp->func))((Pointer)&ompc);
			switch (ret)
			{
				case OM_MSG_RET_ERROR:
					ErrShow();
					break;
				case OM_MSG_RET_DEL:
					break;
				case OM_MSG_RET_OK:
					break;
				case OM_MSG_RET_DONE:
					sip = (SeqIdPtr)(ompc.output_data);
					if (sip != NULL)
					{
						if (sip->choice == SEQID_GI)
						{
							gi = sip->data.intvalue;
							SeqIdFree(sip);
							return gi;
						}
						SeqIdFree(sip);
					}
					break;
				default:
					break;
			}
		}
	}

	return gi;
}


/*****************************************************************************
*
*   GetSeqIdForGI(Int4)
*     returns the SeqId for a GI
*     returns NULL if can't find it
*     The returned SeqId is allocated. Caller must free it.
*
*****************************************************************************/
typedef struct seqidblock {
  Int4       uid;
  SeqIdPtr   sip;
} SeqIdBlock, PNTR SeqIdBlockPtr;

static ValNodePtr seqidgicache = NULL;
static ValNodePtr PNTR seqidgiarray = NULL;
static Int2 seqidcount = 0;
static TNlmRWlock sid_RWlock = NULL;

static void RecordInSeqIdGiCache (Int4 gi, SeqIdPtr sip)

{
	ValNodePtr vnp, tmp;
	ValNodePtr PNTR prev;
	SeqIdBlockPtr sibp;
	Int4 retval;

	/* if (sip == NULL) return; okay to cache NULL because we protect against SeqIdDup */

	retval = NlmRWwrlock(sgi_RWlock);
	if (retval != 0)
	{
		ErrPostEx(SEV_ERROR,0,0,"RecSeqIdGi: RWwrlock error [%ld]",
			(long)retval);
		return;
	}

	
	vnp = ValNodeNew (NULL);
	if (vnp == NULL) goto ret;
	sibp = (SeqIdBlockPtr) MemNew (sizeof (SeqIdBlock));
	if (sibp == NULL) {
		MemFree (vnp);
		goto ret;
	}

	sibp->uid = gi;
	if (sip != NULL) {
		sibp->sip = SeqIdDup (sip);
	}
	vnp->data.ptrvalue = (Pointer) sibp;

	if (seqidgicache == NULL) {
		seqidgicache = vnp;
		goto ret;
	}

	seqidgiarray = MemFree (seqidgiarray);

	prev = (ValNodePtr PNTR) (&seqidgicache);
	tmp = seqidgicache;
	while (tmp != NULL) {
		sibp = (SeqIdBlockPtr) tmp->data.ptrvalue;
		if (sibp != NULL) {
			if (sibp->uid > gi) {
				if (prev != NULL) {
					vnp->next = *prev;
					*prev = vnp;
				}
				goto ret;
			} else if (sibp->uid == gi) {
				goto ret;
			} else {
				prev = (ValNodePtr PNTR) (& (tmp->next));
			}
		} else {
			prev = (ValNodePtr PNTR) (& (tmp->next));
		}
		tmp = tmp->next;
	}
	if (prev != NULL) {
		vnp->next = *prev;
		*prev = vnp;
	}
ret:
	retval = NlmRWunlock(sgi_RWlock);
	if (retval != 0)
	{
		ErrPostEx(SEV_ERROR,0,0,"RecSeqIdGiUnlock: RWunlock error [%ld]",
			(long)retval);
	}
	return;
}

static Boolean FetchFromSeqIdGiCache (Int4 gi, SeqIdPtr PNTR sipp)

{
	ValNodePtr vnp;
	SeqIdBlockPtr sibp = NULL;
	Int2 i;
	Int2 left, right, mid;
	Int4 compare, ret;
	Boolean done = FALSE;

	if (sipp == NULL) return done;
	*sipp = NULL;
	if (seqidgicache == NULL) return done;

	if (seqidgiarray == NULL) {
	    ret = NlmRWwrlock(sgi_RWlock);
		if (ret != 0)
		{
			ErrPostEx(SEV_ERROR,0,0,"SeqIdGi: RWwrlock error [%ld]",
				(long)ret);
			return done;
		}

		if (seqidgiarray == NULL)
		{
			seqidcount = (Int2) ValNodeLen (seqidgicache);
			seqidgiarray = MemNew (sizeof (ValNodePtr) * (size_t) (seqidcount + 1));
			if (seqidgiarray != NULL) {
				i = 0;
				for (vnp = seqidgicache; vnp != NULL; vnp = vnp->next) {
					seqidgiarray [i] = vnp;
					i++;
				}
			}
		}
		ret = NlmRWunlock(sgi_RWlock);
		if (ret != 0)
		{
			ErrPostEx(SEV_ERROR,0,0,"SeqIdGiUnlock: RWunlock error [%ld]",
				(long)ret);
			return done;
		}

	}

	ret = NlmRWrdlock(sgi_RWlock);
    if (ret != 0)
	{
		ErrPostEx(SEV_ERROR,0,0,"SeqIdGi: RWrdlock error [%ld]",
			(long)ret);
		return done;
	}

	if (seqidgiarray != NULL) {
		left = 1;
		right = seqidcount;
		while (left <= right) {
			mid = (left + right) / 2;
			compare = 0;
			vnp = seqidgiarray [mid - 1];
			if (vnp != NULL) {
				sibp = (SeqIdBlockPtr) vnp->data.ptrvalue;
				if (sibp != NULL) {
					compare = gi - sibp->uid;
				}
			}
			if (compare <= 0) {
				right = mid - 1;
			}
			if (compare >= 0) {
				left = mid + 1;
			}
		}
		if (left > right + 1 && sibp != NULL) {
			if (sibp->sip != NULL) {
				*sipp = SeqIdDup (sibp->sip);
			}
			done = TRUE;
		}
	}


	ret = NlmRWunlock(sgi_RWlock);
	if (ret != 0)
	{
		ErrPostEx(SEV_ERROR,0,0,"SeqIdGiUnlock: RWunlock error [%ld]",
			(long)ret);
	}

	return done;
}

NLM_EXTERN SeqIdPtr LIBCALL GetSeqIdForGI (Int4 gi)
{
	BioseqPtr bsp = NULL;
	ObjMgrProcPtr ompp;
	OMProcControl ompc;
	Int2 ret;
	SeqIdPtr sip, sip2=NULL, other=NULL, gb=NULL;
	ValNode vn;


	if (gi <= 0)
		return sip2;

	vn.choice = SEQID_GI;
	vn.data.intvalue = gi;
	vn.next = NULL;

	bsp = BioseqFindCore(&vn);
	if (bsp != NULL)
	{
		for (sip = bsp->id; sip != NULL; sip = sip->next)
		{
		    switch (sip->choice) 
		    {
		        case SEQID_LOCAL:           /* object id */
      		  case SEQID_GIBBSQ:         
		        case SEQID_GIBBMT:
      		  case SEQID_PATENT:
		        case SEQID_GENERAL:
						other = sip;
						break;
					case SEQID_GI:
					   break;
		        case SEQID_GENBANK:
      		  case SEQID_EMBL:
		        case SEQID_PIR:
      		  case SEQID_SWISSPROT:
		        case SEQID_DDBJ:
				case SEQID_PRF:
				  case SEQID_PDB:
						gb = sip;
						break;
					default:
						if (other == NULL)
							other = sip;
						break;
		    }
		}
	}


	if (gb != NULL)
		sip2 = gb;
	else if (other != NULL)
		sip2 = other;

	if (sip2 != NULL)
		return SeqIdDup(sip2);

	if (FetchFromSeqIdGiCache (gi, &sip2)) {
		return sip2;
	}

	ompp = NULL;
	while ((ompp = ObjMgrProcFindNext(NULL, OMPROC_FETCH, OBJ_SEQID, OBJ_SEQID, ompp)) != NULL)
	{
		if ((ompp->subinputtype == SEQID_GI) && (ompp->suboutputtype == 0))
		{
			MemSet((Pointer)(&ompc), 0, sizeof(OMProcControl));
			ompc.input_data = &vn;
			ompc.proc = ompp;
			ret = (* (ompp->func))((Pointer)&ompc);
			switch (ret)
			{
				case OM_MSG_RET_ERROR:
					ErrShow();
					break;
				case OM_MSG_RET_DEL:
					break;
				case OM_MSG_RET_OK:
					break;
				case OM_MSG_RET_DONE:
					sip2 = (SeqIdPtr)(ompc.output_data);
					if (sip2 != NULL)
                                          {
                                            RecordInSeqIdGiCache (gi, sip2);
                                            return sip2;
                                          }
					break;
				default:
					break;
			}
		}
	}

	RecordInSeqIdGiCache (gi, sip2);
	return sip2;
}

/*****************************************************************************
*
*   SeqEntryFind(sip)
*   	returns top level seqentry for sip
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr LIBCALL SeqEntryFind (SeqIdPtr sid)
{
	BioseqPtr bsp;
	ObjMgrDataPtr omdp;
	ObjMgrDataPtr PNTR omdpp;
	SeqEntryPtr result=NULL;
	SeqSubmitPtr ssp;
	Int2 i;
	ObjMgrPtr omp;

	bsp = BioseqFind(sid);
	if (bsp == NULL) return result;

	omp = ObjMgrReadLock();
	omdpp = omp->datalist;

	i = ObjMgrLookup(omp, (Pointer)bsp);
	omdp = omdpp[i];
	while (omdp->parentptr != NULL)
	{
		i = ObjMgrLookup(omp, (omdp->parentptr));
		omdp = omdpp[i];
	}

	if (omdp->datatype == OBJ_SEQSUB) {
		ssp = (SeqSubmitPtr) omdp->dataptr;
		if (ssp != NULL && ssp->datatype == 1) {
			result = (SeqEntryPtr) ssp->data;
		}
	} else {
		result = omdp->choice;
	}
	ObjMgrUnlock();
	return result;
}

/*****************************************************************************
*
*   BioseqContextPtr BioseqContextNew (bsp)
*
*****************************************************************************/
NLM_EXTERN BioseqContextPtr LIBCALL BioseqContextNew (BioseqPtr bsp)
{
	ObjMgrDataPtr omdp;
	ObjMgrDataPtr PNTR omdpp;
	Int2 i, ctr=0;
	SeqEntryPtr seps[BIOSEQCONTEXTMAX];
	BioseqContextPtr bcp;
	ObjMgrPtr omp;

	if (bsp == NULL)
		return NULL;


	bcp = MemNew(sizeof(BioseqContext));
	bcp->bsp = bsp;
	bcp->se.choice = 1;   /* bioseq */
	bcp->se.data.ptrvalue = bsp;

	omp = ObjMgrReadLock();
	if (omp == NULL) return BioseqContextFree(bcp);
	omdpp = omp->datalist;

	i = ObjMgrLookup(omp, (Pointer)bsp);
	omdp = omdpp[i];

	if (omdp->choice != NULL)
	{
		seps[ctr] = omdp->choice;
		ctr++;

		while (omdp->parentptr != NULL)
		{
			i = ObjMgrLookup(omp, (omdp->parentptr));
			omdp = omdpp[i];
			if (omdp->choice != NULL)
			{
				if (ctr == BIOSEQCONTEXTMAX)
					ErrPostEx(SEV_ERROR, 0,0, "BioseqContextNew: more than %d levels",(int)ctr);
				else
				{
					seps[ctr] = omdp->choice;
					ctr++;
			    }
			}
		}

		bcp->count = ctr;
		for (i = 0; i < bcp->count; i++)
		{
			ctr--;
			bcp->context[i] = seps[ctr];
		}
	}

	if (omdp->tempload == TL_CACHED)
	{
		ErrPostEx(SEV_ERROR,0,0,"BioseqContextNew: bsp is TL_CACHED");
		bcp = BioseqContextFree(bcp);
	}

	ObjMgrUnlock();

	return bcp;
}

/*****************************************************************************
*
*   BioseqContextFree(bcp)
*
*****************************************************************************/
NLM_EXTERN BioseqContextPtr LIBCALL BioseqContextFree(BioseqContextPtr bcp)
{
	return MemFree(bcp);
}

/*****************************************************************************
*
*   BioseqContextGetSeqDescr(bcp, type, curr, SeqEntryPtr PNTR sep)
*       returns pointer to the next SeqDescr of this type
*       type gives type of Seq-descr
*       if (type == 0)
*          get them all
*       curr is NULL or previous node of this type found
*       moves up from bsp
*		if (sep != NULL) sep set to SeqEntryPtr containing the SeqDescr.
*
*****************************************************************************/
NLM_EXTERN ValNodePtr LIBCALL BioseqContextGetSeqDescr (BioseqContextPtr bcp, Int2 type, ValNodePtr curr, SeqEntryPtr PNTR the_sep)    /* the last one you used */
{
	Int2 i;
	ValNodePtr tmp;
	Boolean found = FALSE;
	BioseqPtr bsp;
	BioseqSetPtr bssp;

	if (bcp == NULL) return NULL;
	
	if (the_sep != NULL)
		*the_sep = NULL;
		
	if (bcp->count == 0)   /* just a Bioseq */
	{
		tmp = BioseqGetSeqDescr(bcp->bsp, type, curr);
		if (the_sep != NULL) *the_sep = bcp->context[1];
		return tmp;
	}

	i = bcp->count - 1;
	if (curr != NULL)   /* find where we are */
	{
		while ((i >= 0) && (! found))
		{
			if (IS_Bioseq(bcp->context[i]))
			{
				bsp = (BioseqPtr)((bcp->context[i])->data.ptrvalue);
				tmp = bsp->descr;
			}
			else
			{
				bssp = (BioseqSetPtr)((bcp->context[i])->data.ptrvalue);
				tmp = bssp->descr;
			}
			while ((tmp != curr) && (tmp != NULL))
				tmp = tmp->next;
			if (tmp == curr)
			{
				found = TRUE;
				tmp = tmp->next;
			}
			else
				i--;
		}
		if (! found)   /* can't find it! */
			return NULL;
	}
	else              /* get first one */
	{
		tmp = bcp->bsp->descr;
	}
		
	while (i >= 0)
	{
		while (tmp != NULL)
		{
			if ((! type) || ((Int2)(tmp->choice) == type))
			{
				if (the_sep != NULL) *the_sep = bcp->context[i];
				return tmp;
			}
			tmp = tmp->next;
		}
		i--;
		if (i >= 0)
		{
			if (IS_Bioseq(bcp->context[i]))
			{
				bsp = (BioseqPtr)((bcp->context[i])->data.ptrvalue);
				tmp = bsp->descr;
			}
			else
			{
				bssp = (BioseqSetPtr)((bcp->context[i])->data.ptrvalue);
				tmp = bssp->descr;
			}
		}
	}
    return NULL;
}

/*****************************************************************************
*
*   BioseqContextGetSeqFeat(bcp, type, curr, sapp)
*       returns pointer to the next Seq-feat of this type
*       type gives type of Seq-descr
*       if (type == 0)
*          get them all
*       curr is NULL or previous node of this type found
*       moves up from bsp
*   	if (sapp != NULL) is filled with SeqAnnotPtr containing the SeqFeat
*   	in:
*   		0 = sfp->location only
*   		1 = sfp->product only
*   		2 = either of above
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL BioseqContextGetSeqFeat (BioseqContextPtr bcp, Int2 type,
	SeqFeatPtr curr, SeqAnnotPtr PNTR sapp, Int2 in)    /* the last one you used */
{
	SeqEntryPtr sep;

	if (bcp == NULL) return NULL;
	
	if (sapp != NULL)
		*sapp = NULL;

	if (bcp->count == 0)    /* just a BioseqSeq */
		sep = &(bcp->se);
	else
		sep = bcp->context[0];

	return SeqEntryGetSeqFeat (sep, type, curr, sapp, in, bcp->bsp);
}

typedef struct smgetseqfeat {
	Boolean hit;
	SeqFeatPtr last,
		this;
	SeqAnnotPtr sap;
	SeqLocPtr slp1, slp2;
	Int2 in, type;
} SMGetSeqFeat, PNTR GetSeqFeatPtr;

NLM_EXTERN void GetSeqFeatCallback (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);

/*****************************************************************************
*
*   SeqEntryGetSeqFeat(sep, type, curr, sapp)
*       returns pointer to the next Seq-feat of this type
*       type gives type of SeqFeat
*       if (type == 0)
*          get them all
*       curr is NULL or previous node of this type found
*       moves up from bsp
*   	if (sapp != NULL) is filled with SeqAnnotPtr containing the SeqFeat
*       if (bsp != NULL) then for that Bioseq match on location by "in"
*   	in:
*   		0 = sfp->location only
*   		1 = sfp->product only
*   		2 = either of above
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL SeqEntryGetSeqFeat (SeqEntryPtr sep, Int2 type,
	SeqFeatPtr curr, SeqAnnotPtr PNTR sapp, Int2 in, BioseqPtr bsp)    /* the last one you used */
{
	SMGetSeqFeat gsf;
	ValNode vn1, vn2;
	
	if (sep == NULL) return NULL;
	
	if (sapp != NULL)
		*sapp = NULL;

	if (curr == NULL)
		gsf.hit = TRUE;
	else
		gsf.hit = FALSE;
	gsf.last = curr;
	gsf.this = NULL;
	gsf.sap = NULL;
	gsf.type = type;
	gsf.in = in;
	if (bsp != NULL)   /* matching by Bioseq */
	{
		if ((bsp->repr == Seq_repr_seg) || (bsp->repr == Seq_repr_ref))
		{
			vn2.choice = SEQLOC_MIX;
			vn2.data.ptrvalue = bsp->seq_ext;
			gsf.slp2 = (SeqLocPtr)(&vn2);
		}
		else
			gsf.slp2 = NULL;

		vn1.choice = SEQLOC_WHOLE;
		vn1.data.ptrvalue = (Pointer) SeqIdFindBest (bsp->id, 0);
		gsf.slp1 = (SeqLocPtr)(&vn1);
	}
	else
		gsf.slp1 = NULL;

	SeqEntryExplore (sep, (Pointer)(&gsf), GetSeqFeatCallback);

	if (sapp != NULL)
		*sapp = gsf.sap;

	return gsf.this;
}

NLM_EXTERN void GetSeqFeatCallback (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	GetSeqFeatPtr gsfp;
	BioseqPtr bsp;
	BioseqSetPtr bssp;
	SeqAnnotPtr sap;
	SeqFeatPtr sfp, last;
	Boolean hit, gotit = FALSE;
	Uint1 type;
	SeqLocPtr slp1, slp2, slp;
	Int2 i, in, retval;

	gsfp = (GetSeqFeatPtr)data;
	if (gsfp->this != NULL)   /* got it */
		return;

	last = gsfp->last;
	hit = gsfp->hit;
	type = (Uint1)(gsfp->type);
	if (gsfp->slp1 != NULL)   /* matching by Bioseq */
	{
		slp1 = gsfp->slp1;
		slp2 = gsfp->slp2;
		in = gsfp->in;
	}
	else
		slp1 = NULL;

	if (IS_Bioseq(sep))
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		sap = bsp->annot;
	}
	else
	{
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		sap = bssp->annot;
	}

	while (sap != NULL)
	{
		if (sap->type == 1)  /* feature table */
		{
			for (sfp = (SeqFeatPtr)(sap->data); sfp != NULL; sfp = sfp->next)
			{
				if (! hit)	   /* still looking */
				{
					if (sfp == last)
					{
						hit = TRUE;
						gsfp->hit = TRUE;
					}
				}
				else
				{
					if ((! type) || (type == sfp->data.choice))
					{
						if (slp1 != NULL)   /* look for feats on a bioseq */
						{
							for (i = 0; i < 2; i++)
							{
								if ((i == 0) && (in != 1))
									slp = sfp->location;
								else if ((i==1) && (in != 0))
									slp = sfp->product;
								else
									slp = NULL;
								if (slp != NULL)
								{
									retval = SeqLocCompare(slp, slp1);
									if (retval)
									{
										gotit = TRUE;
										break;
									}

									if (slp2 != NULL)
									{
										retval = SeqLocCompare(slp, slp2);
										if (retval)
										{
											gotit = TRUE;
											break;
										}
									}
								}
							}
						}
						else
							gotit = TRUE;
						if (gotit)
						{
							gsfp->this = sfp;
							gsfp->sap = sap;
							return;
						}
					}
				}
			}
		}

		sap = sap->next;
	}

	return;
}

/*****************************************************************************
*
*   BioseqContextGetTitle(bcp)
*     returns first title for Bioseq in this context
*
*****************************************************************************/
NLM_EXTERN CharPtr LIBCALL BioseqContextGetTitle(BioseqContextPtr bcp)
{
	CharPtr title = NULL;
	ValNodePtr vnp;

	vnp = BioseqContextGetSeqDescr(bcp, Seq_descr_title, NULL, NULL);
	if (vnp != NULL)
		title = (CharPtr)vnp->data.ptrvalue;
	return title;
}

/*****************************************************************************
*
*   SeqMgr Functions
*
*****************************************************************************/

/*****************************************************************************
*
*   SeqMgrSeqEntry(type, data, sep)
*   	Adds the SeqEntryPtr pointing directly to this Bioseq or BioseqSet
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrSeqEntry (Uint1 type, Pointer data, SeqEntryPtr sep)
{
	return ObjMgrSetChoice (OBJ_SEQENTRY, sep, data);
}

/*****************************************************************************
*
*   SeqMgrGetSeqEntryForData(data)
*   	returns SeqEntryPtr for a BioseqPtr or BioseqSetPtr
*       sep must have been put in SeqMgr using SeqMgrSeqEntry
*       the Bioseq/BioseqSets it is a part of must also be in SeqMgr
*       returns NULL on failure.
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr LIBCALL SeqMgrGetSeqEntryForData (Pointer data)
{
	return ObjMgrGetChoiceForData(data);
}

/*****************************************************************************
*
*   SeqMgrGetEntityIDForSeqEntry(sep)
*   	returns the EntityID for a SeqEntryPtr
*       sep must have been put in SeqMgr using SeqMgrSeqEntry
*       the Bioseq/BioseqSets it is a part of must also be in SeqMgr
*       This function will move up to the top of the SeqEntry tree it may be
*          in. If top level EntityID is 0, one is assigned at this point.
*       If an element is moved under a different hierarchy, its EntityID will
*          change.
*       returns 0 on failure.
*
*****************************************************************************/
NLM_EXTERN Int2 LIBCALL SeqMgrGetEntityIDForSeqEntry (SeqEntryPtr sep)
{
	return ObjMgrGetEntityIDForChoice (sep);
}

/*****************************************************************************
*
*   SeqMgrGetSeqEntryForEntityID (id)
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr LIBCALL SeqMgrGetSeqEntryForEntityID (Int2 id)
{
	return ObjMgrGetChoiceForEntityID (id);
}

/*****************************************************************************
*
*   SeqMgrSetBSFetchTop (fetch, data)
*   	sets the BSFetchTop routine to "fetch"
*       returns previous value
*       set to NULL to turn off all fetching for that type
*
*       Current value can be called directly as BioseqFetch();
*   	Default is
*   		1) looks in memory
*   		2) looks locally if LocalBSFetch is set
*			3) looks remotely if RemoteBSFetch is set
*
*****************************************************************************/
NLM_EXTERN BSFetchTop LIBCALL SeqMgrSetBSFetchTop (BSFetchTop fetch, Pointer data)
{
	SeqMgrPtr smp;
	BSFetchTop tmp = NULL;

	smp = SeqMgrWriteLock();
	if (smp == NULL) return tmp;
	
	tmp = smp->bsfetch;
	smp->bsfetch = fetch;
	smp->TopData = data;
	SeqMgrUnlock();
	return tmp;
}

/*****************************************************************************
*
*   SeqMgrSetFetchOnLock(value)
*   	if value = TRUE, manager will try to fetch the bioseq if not in
*          memory, when BioseqLock is called
*   	if FALSE, BioseqLock will only look in memory
*       returns previous value of fetch_on_lock
*       default is to fetch_on_lock
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrSetFetchOnLock (Boolean value)
{
	Boolean tmp=FALSE;
	SeqMgrPtr smp;

	smp = SeqMgrWriteLock();
	if (smp == NULL) return tmp;

	tmp = smp->fetch_on_lock;
	smp->fetch_on_lock = value;
	SeqMgrUnlock();
	return tmp;
}

/*****************************************************************************
*
*   SeqMgrLinkSeqEntry(sep, parenttype, parentptr)
*      connects all component seq-entries by traversing the linked list
*        all calling SeqMgrConnect and SeqMgrSeqEntry appropriately
*        if parenttype != 0, then assumes seqentry is contained in parentptr
*           and should be connected to it
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrLinkSeqEntry (SeqEntryPtr sep, Uint2 parenttype, Pointer parentptr)
{
	SeqEntryPtr sep2;
	BioseqSetPtr bssp;
	Uint2 the_type;
	
	if (sep == NULL)
		return FALSE;

	if (IS_Bioseq(sep))
		the_type = OBJ_BIOSEQ;
	else
		the_type = OBJ_BIOSEQSET;

	SeqMgrSeqEntry((Uint1)the_type, sep->data.ptrvalue, sep);

	/**** if (parenttype != 0) ****/
	ObjMgrConnect(the_type, sep->data.ptrvalue, parenttype, parentptr);

	if (! IS_Bioseq(sep))
	{
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		for (sep2 = bssp->seq_set; sep2 != NULL; sep2 = sep2->next)
		{
			SeqMgrLinkSeqEntry (sep2, the_type, sep->data.ptrvalue);
		}
	}
	return TRUE;
}
/*****************************************************************************
*
*   Selection Functions for data objects based on SeqLoc
*      See also general selection in objmgr.h
*
*****************************************************************************/

/*****************************************************************************
*
*   SeqMgrSelect(region)
*      region is a SeqLocPtr
*          It can only apply to one Bioseq
*          selected area will be extreme left and right ends
*          fuzziness is ignored
*      if something else selected, deselects it first, then selects requested
*        item
*      to select without deselecting something else, use SeqMgrAlsoSelect()
*      returns TRUE if item is now currently selected, FALSE if not
*      "region" is always copied. Caller is responsible for managment of
*         SeqLoc that is passed in.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrSelect (SeqLocPtr region)
{
	return SeqMgrGenericSelect(region, 1, NULL);
}

NLM_EXTERN Boolean LIBCALL SeqMgrAlsoSelect (SeqLocPtr region)
{
	return SeqMgrGenericSelect(region, 2, NULL);
}

/*****************************************************************************
*
*   SeqMgrDeSelect(region)
*   	if this item was selected, then deselects and returns TRUE
*   	else returns FALSE
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrDeSelect (SeqLocPtr region)
{
	return SeqMgrGenericSelect(region, 3, NULL);
}

/*****************************************************************************
*
*   SeqMgrSetColor(region, rgb)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrSetColor (SeqLocPtr region, Uint1Ptr rgb)
{
	if (rgb == NULL) return FALSE;
        return SeqMgrGenericSelect(region, 4, rgb);
}

static Boolean NEAR SeqMgrGenericSelect (SeqLocPtr region, Int2 type,
                                           Uint1Ptr rgb)
{
	SeqInt si;
	ValNode vn;
	SeqIdPtr sip;
	Uint2 entityID, itemID;

	if (region == NULL) return FALSE;

	sip = SeqLocId(region);
	if (sip == NULL) return FALSE;

	entityID = BioseqFindEntity(sip, &itemID);
	if (entityID == 0) return FALSE;

	MemSet((Pointer)(&si), 0, sizeof(SeqInt));
	MemSet((Pointer)(&vn), 0, sizeof(ValNode));

	si.id = sip;
	si.from = SeqLocStart(region);
	si.to = SeqLocStop(region);
	si.strand = SeqLocStrand(region);

	if ((si.from < 0) || (si.to < 0) || (si.from > si.to)) return FALSE;

	vn.choice = SEQLOC_INT;
	vn.data.ptrvalue = &si;

	switch (type)
	{
		case 1:
			return ObjMgrSelect(entityID, itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, &vn);
		case 2:
			return ObjMgrAlsoSelect(entityID, itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, &vn);
		case 3:
			return ObjMgrDeSelect(entityID, itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, &vn);
		case 4:
			return ObjMgrSetColor(entityID, itemID, OBJ_BIOSEQ,
                                 OM_REGION_SEQLOC, &vn, rgb);
		default:
			break;
	}

	return FALSE;
}

/*****************************************************************************
*
*   SpreadGapsInDeltaSeq(BioseqPtr bsp)
*      bsp must be a delta seq
*      function counts deltas with known lengths ( = known_len)
*               counts deltas which are gaps of unknown length ( = unk_count)
*                  these can delta of length 0, delta with fuzz = lim (unk),
*                    or SEQLOC_NULL
*               converts all unknown gaps to delta with fuzz = lim(unk)
*               sets length of all unknown gaps to
*                  (bsp->length - known_len)/unk_count
*                  any reminder spread over first few gaps
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SpreadGapsInDeltaSeq (BioseqPtr bsp)
{
	Boolean retval = FALSE;
	Int4 known_len = 0,
		 total_gap, gap_len,
		 unk_count = 0,
		 remainder;
	DeltaSeqPtr dsp;
	SeqLocPtr slocp;
	SeqLitPtr slp;
	IntFuzzPtr ifp;

	if (bsp == NULL) return retval;
	if ((bsp->repr != Seq_repr_delta) || (bsp->seq_ext == NULL))
		return retval;

	retval = TRUE;  /* can function */

	for (dsp = (DeltaSeqPtr)(bsp->seq_ext); dsp != NULL; dsp = dsp->next)
	{
		switch (dsp->choice)
		{
			case 1:	  /* SeqLocPtr */
				slocp = (SeqLocPtr)(dsp->data.ptrvalue);
				if (slocp == NULL) break;
				if (slocp->choice == SEQLOC_NULL)  /* convert it */
				{
					SeqLocFree(slocp);
					slp = SeqLitNew();
					dsp->choice = 2;
					dsp->data.ptrvalue = slp;
					ifp = IntFuzzNew();
					slp->fuzz = ifp;
					ifp->choice = 4;   /* lim - type unk */
					unk_count++;
				}
				else                               /* count length */
					known_len += SeqLocLen(slocp);
				break;
			case 2:   /* SeqLitPtr */
				slp = (SeqLitPtr)(dsp->data.ptrvalue);
				if (slp == NULL) break;
				if (slp->seq_data != NULL)         /* not a gap */
				{
					known_len += slp->length;
					break;
				}
				ifp = slp->fuzz;
				if (slp->length == 0)  /* unknown length */
				{
					unk_count++;
					if (ifp != NULL)
					{
						if (ifp->choice != 4)  /* not lim */
							ifp = IntFuzzFree(ifp);
						else if (ifp->a != 0)  /* not unk */
							ifp = IntFuzzFree(ifp);
					}
					if (ifp == NULL)
					{
						ifp = IntFuzzNew();
						ifp->choice = 4; /* lim - unk */
						slp->fuzz = ifp;
					}
				}
				else                      /* gap length was set */
				{
					if (ifp == NULL)  /* no fuzz - count length */
						known_len += slp->length;
					else              /* might be a guess */
					{
						if ((ifp->choice == 4) && (ifp->a == 0)) /* lim - unk */
							unk_count++;
						else
							known_len += slp->length;
					}
				}
				break;
			default:
				break;
		}

	}

	if (unk_count == 0)   /* no unknown gaps */
		return retval;

	total_gap = bsp->length - known_len;
	if (total_gap < 0)
		total_gap = 0;
	gap_len = total_gap / unk_count;
	remainder = total_gap - (gap_len * unk_count);

	for (dsp = (DeltaSeqPtr)(bsp->seq_ext); dsp != NULL; dsp = dsp->next)
	{
		switch (dsp->choice)
		{
			case 1:	  /* SeqLocPtr */
				break;
			case 2:   /* SeqLitPtr */
				slp = (SeqLitPtr)(dsp->data.ptrvalue);
				if (slp == NULL) break;
				if (slp->seq_data != NULL) break;
				ifp = slp->fuzz;
				if (ifp == NULL) break;
				if ((ifp->choice != 4) || (ifp->a != 0))
					break;
				slp->length = gap_len;
				if (remainder)
				{
					slp->length++;
					remainder--;
				}
				break;
			default:
				break;
		}
	}

	return retval;
}

/*****************************************************************************
*
*   CountGapsInDeltaSeq(BioseqPtr bsp, &num_segs, &num_gaps, &known_residues, &num_gaps_faked)
*      bsp must be a delta seq
*      function counts deltas and returns a profile
*          num_segs = total number of segments
*          num_gaps = total number of segments representing gaps
*          known_residues = number of real residues in the sequence (not gaps)
*          num_gaps_faked = number of gaps where real length is not known, but where
*                           a length was guessed by spreading the total gap length
*                           out over all gaps evenly.
*
*      NOTE: any of these pointers except bsp can be NULL
*
*      returns TRUE if values in argument were set.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL CountGapsInDeltaSeq (BioseqPtr bsp, Int4Ptr num_segs, Int4Ptr num_gaps, Int4Ptr known_residues, Int4Ptr num_gaps_faked, CharPtr buf, Int2 buflen)
{
	Boolean retval = FALSE;
	Int4 residues = 0,
		segs = 0,
		gaps = 0,
		len = 0,
		fake_gaps = 0,
		from = 0, 
		tlen = 0;
	DeltaSeqPtr dsp;
	SeqLocPtr slocp;
	SeqLitPtr slp;
	IntFuzzPtr ifp;
	Boolean unk;
	static Char tmp[128];
	Int2 diff;

	if (bsp == NULL) return retval;
	if ((bsp->repr != Seq_repr_delta) || (bsp->seq_ext == NULL))
		return retval;

	retval = TRUE;  /* can function */


	for (dsp = (DeltaSeqPtr)(bsp->seq_ext); dsp != NULL; dsp = dsp->next)
	{
		segs++;
		from = len + 1;
		switch (dsp->choice)
		{
			case 1:	  /* SeqLocPtr */
				slocp = (SeqLocPtr)(dsp->data.ptrvalue);
				if (slocp == NULL) break;
				if (slocp->choice == SEQLOC_NULL)  /* gap */
				{
					gaps++;
					sprintf(tmp, "* %ld %ld gap of unknown length~", from, len);
					diff = LabelCopy(buf, tmp, buflen);
					buflen -= diff;
					buf += diff;
				}
				else {                              /* count length */
					residues += SeqLocLen(slocp);
					if (buf != NULL) {
						tlen =  SeqLocLen(slocp);
						len  += tlen;
						sprintf(tmp, "* %8ld %8ld: contig of %ld bp in length~", from, len, tlen);
						diff = LabelCopy(buf, tmp, buflen);
						buflen -= diff;
						buf += diff;
					}
				}
				break;
			case 2:   /* SeqLitPtr */
				slp = (SeqLitPtr)(dsp->data.ptrvalue);
				if (slp == NULL) break;
				tlen =  slp->length;
				len  += tlen;
				if (slp->seq_data != NULL)
				{
					residues += slp->length;
					if (buf) {
						sprintf(tmp, "* %8ld %8ld: contig of %ld bp in length~", from, len, tlen);
						diff = LabelCopy(buf, tmp, buflen);
						buflen -= diff;
						buf += diff;
					}
				}
				else
				{
					unk = FALSE;
					gaps++;
					ifp = slp->fuzz;
					if (ifp != NULL)
					{
						if ((ifp->choice == 4) && (ifp->a == 0)) {
							unk = TRUE;
							fake_gaps++;
							if (buf) {
								if (from > len) {
								sprintf(tmp, "*                    gap of unknown length~");
								} else {
								sprintf(tmp, "* %8ld %8ld: gap of unknown length~", from, len);
								}
								diff = LabelCopy(buf, tmp, buflen);
								buflen -= diff;
								buf += diff;
							}
						}
					}
					if (!unk && buf) {
						sprintf(tmp, "* %8ld %ld: gap of %8ld bp~", from, len, tlen);
						diff = LabelCopy(buf, tmp, buflen);
						buflen -= diff;
						buf += diff;
					}
				}
				break;
			default:
				break;
		}
	}

	if (num_segs != NULL)
		*num_segs = segs;
	if (num_gaps != NULL)
		*num_gaps = gaps;
	if (known_residues != NULL)
		*known_residues = residues;
	if (num_gaps_faked != NULL)
		*num_gaps_faked = fake_gaps;

	return retval;
}


/*****************************************************************************
*
*   SeqMgrAdd(type, data)
*   	adds a Bioseq or BioseqSet to the sequence manager
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrAdd (Uint2 type, Pointer data)
{
	Boolean retval;

	retval = ObjMgrAdd(type, data);
	if (type != OBJ_BIOSEQ)
		return retval;

	SeqMgrAddToBioseqIndex((BioseqPtr)data);

	return retval;

}


/*****************************************************************************
*
*   SeqMgrDelete(type, data)
*   	deletes a Bioseq or BioseqSet from the sequence manager
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrDelete (Uint2 type, Pointer data)
{
	if (type == OBJ_BIOSEQ)  /* remove id indexes */
		SeqMgrDeleteFromBioseqIndex((BioseqPtr)data);

	return ObjMgrDelete(type, data);
}

static Boolean NEAR SeqMgrAddIndexElement(SeqMgrPtr smp, BioseqPtr bsp, CharPtr buf)
{
	SeqIdIndexElementPtr sip, PNTR sipp;
	SeqIdIndexBlockPtr sibp, prev;
	Int4 imin, imax, i, j;
	CharPtr tmp, newstr;
	ObjMgrDataPtr omdp;
	ObjMgrPtr omp;

	omp = ObjMgrReadLock();
	omdp = ObjMgrFindByData(omp, (Pointer)bsp);  /* caching protection */
	ObjMgrUnlock();
	if (omdp == NULL)
	{
		return FALSE;
	}

	sipp = smp->BioseqIndex;
	if (smp->BioseqIndexCnt >= smp->BioseqIndexNum)  /* expand space */
	{
	   prev = NULL;
	   for (sibp = smp->BioseqIndexData; sibp != NULL; sibp = sibp->next)
		   prev = sibp;
	   sibp = MemNew(sizeof(SeqIdIndexBlock));
	   if (prev != NULL)
		   prev->next = sibp;
	   else
		   smp->BioseqIndexData = sibp;

	   smp->BioseqIndex = MemNew((smp->BioseqIndexNum + 100) * 
sizeof(SeqIdIndexElementPtr));
	   MemCopy(smp->BioseqIndex, sipp, (smp->BioseqIndexNum * 
sizeof(SeqIdIndexElementPtr)));
	   MemFree(sipp);
	   smp->BioseqIndexNum += 100;
	   sipp = smp->BioseqIndex;
	   for (i = 0, j = smp->BioseqIndexCnt; i < 100; i++, j++)
		   sipp[j] = &(sibp->sid[i]);
	}

	i = smp->BioseqIndexCnt;   /* empties are at the end */
	sip = sipp[i];
	sip->omdp = omdp;       /* fill in the values */
	sip->str = StringSave(buf);
	newstr = sip->str;
	RevStringUpper(newstr);  /* try to avoid case check */

	imin = 0;                   /* find where it goes */
	imax = i-1;
	if (imax >= 0)
		tmp = sipp[imax]->str;
	if ((i) && (StringCmp(newstr, sipp[imax]->str) < 0))
	{
		i = (imax + imin) / 2;
		while (imax > imin)
		{
			tmp = sipp[i]->str;
			if ((j = StringCmp(newstr, tmp)) < 0)
				imax = i - 1;
			else if (j > 0)
				imin = i + 1;
			else
				break;
			i = (imax + imin)/2;
		}

		if (StringCmp(newstr, sipp[i]->str) > 0) /* check for off by 1 */
		{
			i++;
		}


		imax = smp->BioseqIndexCnt - 1;	 /* open the array */
		while (imax >= i)
		{
			sipp[imax+1] = sipp[imax];
			imax--;
		}
	}

	sipp[i] = sip;    /* put in the pointer in order */
	smp->BioseqIndexCnt++;     /* got one more */
	return TRUE;
}

/*****************************************************************************
*
*   SeqMgrProcessNonIndexedBioseq()
*   	Indexes a BioseqPtr by SeqId(s)
*
*****************************************************************************/
static Boolean NEAR SeqMgrProcessNonIndexedBioseq(void)
{
	BioseqPtr PNTR bspp, bsp;
	Int4 i, total, k;
	SeqIdPtr sip;
	Char buf[80];
	CharPtr tmp;
	Uint1 oldchoice;
	Boolean indexed;
	TextSeqIdPtr tsip;
	SeqMgrPtr smp;

	smp = SeqMgrReadLock();
	if (! smp->NonIndexedBioseqCnt)
	{
		SeqMgrUnlock();
		return TRUE;
	}
	SeqMgrUnlock();

	smp = SeqMgrWriteLock();
	if (! smp->NonIndexedBioseqCnt)
	{
		SeqMgrUnlock();
		return TRUE;
	}

	total = smp->NonIndexedBioseqCnt;
	bspp = smp->NonIndexedBioseq;
	for (i = 0; i < total; i++)
	{
		indexed = FALSE;
		bsp = bspp[i];
		if (bsp != NULL)
		{
			if (bsp->id != NULL)
			{
				indexed = TRUE;
				for (sip = bsp->id; sip != NULL; sip = sip->next)
				{
					oldchoice = 0;
					switch (sip->choice)
					{
					case SEQID_GI:
						sprintf(buf, "%ld", (long)(sip->data.ptrvalue));
						SeqMgrAddIndexElement(smp, bsp, buf);
						break;
					case SEQID_EMBL:
					case SEQID_DDBJ:
						oldchoice = sip->choice;
						sip->choice = SEQID_GENBANK;
					case SEQID_GENBANK:
					case SEQID_PIR:
					case SEQID_OTHER:
					case SEQID_SWISSPROT:
					case SEQID_PRF:
						tsip = (TextSeqIdPtr)(sip->data.ptrvalue);
						if (tsip->name != NULL)
						{
							tmp = tsip->accession;
							tsip->accession = NULL;
							SeqIdWrite(sip, buf, PRINTID_FASTA_SHORT, 79);
							SeqMgrAddIndexElement(smp, bsp, buf);
							tsip->accession = tmp;
						}
						tmp = tsip->name;
						tsip->name = NULL;
 						SeqIdWrite(sip, buf, PRINTID_FASTA_SHORT, 79);
						SeqMgrAddIndexElement(smp, bsp, buf);
						tsip->name = tmp;
						if (oldchoice)
							sip->choice = oldchoice;
						break;
					default:
  						SeqIdWrite(sip, buf, PRINTID_FASTA_SHORT, 79);
						SeqMgrAddIndexElement(smp, bsp, buf);
						break;
					}
				}
			}
		}
		if (indexed)
			bspp[i] = NULL;
	}

	for (i = 0; i < total; i++)
	{
		if (bspp[i] == NULL)
		{
		   total--;
		   for (k = i; k < total; k++)
			   bspp[k] = bspp[k+1];
		   i--;
		}
	}

	smp->NonIndexedBioseqCnt = total;

	SeqMgrUnlock();

	return TRUE;
}



/*****************************************************************************
*
*   SeqMgrAddToBioseqIndex(bsp)
*   	Indexes a BioseqPtr by SeqId(s)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrAddToBioseqIndex (BioseqPtr bsp)
{
	SeqMgrPtr smp;
	BioseqPtr PNTR bspp;

	if (bsp == NULL)
		return FALSE;

	smp = SeqMgrWriteLock();
							  /* increase array as 
needed */
	if (smp->NonIndexedBioseqCnt >= smp->NonIndexedBioseqNum)
	{
		bspp = smp->NonIndexedBioseq;
		smp->NonIndexedBioseq = MemNew((smp->NonIndexedBioseqNum + 10) * 
sizeof (BioseqPtr));
		MemCopy(smp->NonIndexedBioseq, bspp, (smp->NonIndexedBioseqNum * 
sizeof(BioseqPtr)));
		MemFree(bspp);
		smp->NonIndexedBioseqNum += 10;
	}

	smp->NonIndexedBioseq[smp->NonIndexedBioseqCnt] = bsp;
	smp->NonIndexedBioseqCnt++;

	SeqMgrUnlock();

	SeqMgrProcessNonIndexedBioseq();

	return TRUE;
}


/*****************************************************************************
*
*   SeqMgrDeleteDeleteFromBioseqIndex(bsp)
*   	Removes index on BioseqPtr SeqIds
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrDeleteFromBioseqIndex (BioseqPtr bsp)
{
	SeqMgrPtr smp;
	SeqIdIndexElementPtr PNTR sipp, sip;
	Int4 i, j, num;
	BioseqPtr PNTR bspp;
	ObjMgrDataPtr omdp;
	ObjMgrPtr omp;

	smp = SeqMgrWriteLock();
								/* check if not 
indexed yet */
	if (smp->NonIndexedBioseqCnt > 0)
	{
		num = smp->NonIndexedBioseqCnt;
		bspp = smp->NonIndexedBioseq;
		for (i = 0; i < num; i++)
		{
			if (bspp[i] == bsp)
			{
				num--;
				for (j = i; j < num; j++)
					 bspp[j] = bspp[j+1];
				i--;
			}
		}
		smp->NonIndexedBioseqCnt = num;
	}

	num = smp->BioseqIndexCnt;
	sipp = smp->BioseqIndex;
	omp = ObjMgrReadLock();
	omdp = ObjMgrFindByData(omp, (Pointer)bsp);
	ObjMgrUnlock();

	for (i = 0; i < BIOSEQ_CACHE_NUM; i++)   /* remove from BioseqFind cache */
	{
		if (omdp_cache[i] == omdp)
		{
			omdp_cache[i] = NULL;
			se_cache[i] = NULL;
		}
	}

	for (i = 0; i < num; i++)
	{
	   if (sipp[i]->omdp == omdp)
	   {
		   sipp[i]->omdp = NULL;
		   sipp[i]->str = MemFree(sipp[i]->str);
		   sip = sipp[i];
		   for (j = i; j < (num-1); j++)
			   sipp[j] = sipp[j+1];
		   sipp[j] = sip;
		   num--; i--;
	   }
	}

	smp->BioseqIndexCnt = num;

	SeqMgrUnlock();

	return TRUE;
}


/*****************************************************************************
*
*   SeqMgrReplaceInBioseqIndex(bsp)
*   	Replaces index on BioseqPtr SeqIds
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrReplaceInBioseqIndex (BioseqPtr bsp)
{
	SeqMgrDeleteFromBioseqIndex(bsp);
	return SeqMgrAddToBioseqIndex(bsp);
}

/*****************************************************************************
*
*   GetUniGeneIDForSeqId(SeqIdPtr)
*     returns the UniGene ID for a SeqId
*     returns 0 if can't find it, or not a legal unigene id
*     This only applies to genomes division of entrez
*
*****************************************************************************/

/*****************************************************************
*
*	IT IS a KLUDGE!! Add 1,000,000 to the unigene id
*
*****************************************************************/
#define KLUDGE_UNIGENE 1000000	/*the kludge offset val add to unigene sequence*/
#define KLUDGE_FlyBase 2000000	/*the kludge offset for FlyBase*/
#define KLUDGE_JACKSON 3000000  /*the kludge offset for the Mouse data*/
#define KLUDGE_JRGP    4000000  /*the kludge offset for the rice data*/
#define KLUDGE_CESC    5000000  /*the kludge offset for the C. elegans data*/
#define KLUDGE_BSNR    6000000  /*the kludge offset for the B. subtilis data*/
#define KLUDGE_HUMGEN  7000000  /*the kludge offset for the Human genomic data*/
#define KLUDGE_YGG     8000000  /*the kludge offset for the yeast data*/
#define KLUDGE_NCBICG  9000000  /*the kludge offset for small genomes*/
#define KLUDGE_MAIZE  10000000  /*the kludge offset for corn*/

NLM_EXTERN Int4 LIBCALL GetUniGeneIDForSeqId (SeqIdPtr sip)
{
	DbtagPtr db_tag;
	ObjectIdPtr oip;

	if (sip == NULL)
		return 0;


	if(sip->choice != SEQID_GENERAL)
		return 0;

	db_tag = sip->data.ptrvalue;
	if(db_tag == NULL || db_tag->db == NULL)
		return 0;

	oip = db_tag->tag;
	if(oip == NULL || oip->id == 0)
		return 0;

	if(StringCmp(db_tag->db, "UNIGENE") == 0)
		return (KLUDGE_UNIGENE+ oip->id);
	if(StringCmp(db_tag->db, "UniGene") == 0)
		return (KLUDGE_UNIGENE+ oip->id);
	if(StringCmp(db_tag->db, "FlyBase") == 0)
		return (KLUDGE_FlyBase+ oip->id);
	if(StringCmp(db_tag->db, "JACKSON") == 0)
		return (KLUDGE_JACKSON+ oip->id);
	if(StringCmp(db_tag->db, "JRGP") == 0)
		return (KLUDGE_JRGP + oip->id);
	if(StringCmp(db_tag->db, "CESC") == 0)
		return (KLUDGE_CESC + oip->id);
	if(StringCmp(db_tag->db, "BSNR") == 0)
		return (KLUDGE_BSNR + oip->id);
        if(StringCmp(db_tag->db, "HUMGEN") == 0)
                return (KLUDGE_HUMGEN + oip->id);
        if(StringCmp(db_tag->db, "YGG") == 0)
                return (KLUDGE_YGG + oip->id);
        if(StringCmp(db_tag->db, "NCBICG") == 0)
                return (KLUDGE_NCBICG + oip->id);
        if(StringCmp(db_tag->db, "MAIZE") == 0)
                return (KLUDGE_MAIZE + oip->id);
	return 0;

}


/*****************************************************************************
*
*   BioseqExtra extensions to preindex for rapid retrieval
*
*****************************************************************************/

/*
*  remaining to be done are mapping tables for rapid coordinate conversion
*  between genome record and parts, genomic DNA and mRNA, and mRNA and protein
*/

static ObjMgrDataPtr SeqMgrGetOmdpForPointer (Pointer ptr)

{
  ObjMgrDataPtr  omdp;
  ObjMgrPtr      omp;

  if (ptr == NULL) return NULL;
  omp = ObjMgrWriteLock ();
  omdp = ObjMgrFindByData (omp, ptr);
  ObjMgrUnlock ();
  return omdp;
}

static SeqEntryPtr SeqMgrGetTopSeqEntryForEntity (Uint2 entityID)

{
  ObjMgrDataPtr  omdp;
  SeqSubmitPtr   ssp;

  omdp = ObjMgrGetData (entityID);
  if (omdp == NULL) return FALSE;
  switch (omdp->datatype) {
    case OBJ_SEQSUB :
      ssp = (SeqSubmitPtr) omdp->dataptr;
      if (ssp != NULL && ssp->datatype == 1) {
        return (SeqEntryPtr) ssp->data;
      }
      break;
    case OBJ_BIOSEQ :
    case OBJ_BIOSEQSET :
      return (SeqEntryPtr) omdp->choice;
    default :
      break;
  }
  return NULL;
}

static Boolean SeqMgrClearBioseqExtraData (ObjMgrDataPtr omdp)

{
  BioseqExtraPtr  bspextra;
  SMFeatBlockPtr  currf;
  SMSeqIdxPtr     currp;
  Int2            i;
  SMFeatItemPtr   itemf;
  SMFeatBlockPtr  nextf;
  SMSeqIdxPtr     nextp;

  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return FALSE;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return FALSE;

  /* free sorted arrays of pointers into data blocks */

  bspextra->featsByID = MemFree (bspextra->featsByID);
  bspextra->featsBySfp = MemFree (bspextra->featsBySfp);
  bspextra->featsByPos = MemFree (bspextra->featsByPos);
  bspextra->genesByPos = MemFree (bspextra->genesByPos);
  bspextra->pubsByPos = MemFree (bspextra->pubsByPos);
  bspextra->orgsByPos = MemFree (bspextra->orgsByPos);

  /* free arrays to speed mapping from parts to segmented bioseq */

  bspextra->partsByLoc = MemFree (bspextra->partsByLoc);
  bspextra->partsBySeqId = MemFree (bspextra->partsBySeqId);

  /* free data blocks of feature information */

  currf = bspextra->featlisthead;
  while (currf != NULL) {
    nextf = currf->next;

    if (currf->data != NULL) {

      /* free allocated label strings within block items */

      for (i = 0; i < currf->index; i++) {
        itemf = &(currf->data [i]);
        MemFree (itemf->label);
        MemFree (itemf->ivals);
      }

      /* free array of SMFeatItems */

      MemFree (currf->data);
    }

    MemFree (currf);
    currf = nextf;
  }

  /* free data blocks of parts to segment mapping information */

  currp = bspextra->segparthead;
  while (currp != NULL) {
    nextp = currp->next;
    SeqLocFree (currp->slp);
    MemFree (currp->seqIdOfPart);
    MemFree (currp);
    currp = nextp;
  }

  /* clean interval list once implemented */

  bspextra->featlisthead = NULL;
  bspextra->featlisttail = NULL;
  bspextra->segparthead = NULL;

  bspextra->numfeats = 0;
  bspextra->numgenes = 0;
  bspextra->numpubs = 0;
  bspextra->numorgs = 0;
  bspextra->numsegs = 0;

  bspextra->min = INT4_MAX;
  bspextra->blocksize = 50;

  bspextra->protFeat = NULL;
  bspextra->cdsOrRnaFeat = NULL;

  /* free genome - parts mapping arrays when they are added */

  return TRUE;
}

static Boolean DoSeqMgrFreeBioseqExtraData (ObjMgrDataPtr omdp)

{
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return FALSE;
  if (omdp->extradata != NULL) {
    SeqMgrClearBioseqExtraData (omdp);
    omdp->extradata = MemFree (omdp->extradata);
    omdp->reapextra = NULL;
    omdp->reloadextra = NULL;
    omdp->freeextra = NULL;
  }
  return TRUE;
}

/* object manager callbacks to reap, reload, and free extra bioseq data */

NLM_EXTERN Pointer LIBCALLBACK SeqMgrReapBioseqExtraFunc (Pointer data)

{
  BioseqExtraPtr  bspextra;
  SMFeatBlockPtr  curr;
  Int2            i;
  SMFeatItemPtr   item;
  ObjMgrDataPtr   omdp;

  omdp = (ObjMgrDataPtr) data;
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;

  /* loop through data blocks of feature information */

  curr = bspextra->featlisthead;
  while (curr != NULL) {

    /* NULL out pointers to cached out feature and annot */

	if (curr->data != NULL) {
      for (i = 0; i < curr->index; i++) {
        item = &(curr->data [i]);
        item->sfp = NULL;
        item->sap = NULL;
      }
    }

    curr = curr->next;
  }

  return NULL;
}

NLM_EXTERN Pointer LIBCALLBACK SeqMgrReloadBioseqExtraFunc (Pointer data)

{
  return NULL;
}

NLM_EXTERN Pointer LIBCALLBACK SeqMgrFreeBioseqExtraFunc (Pointer data)

{
  DoSeqMgrFreeBioseqExtraData ((ObjMgrDataPtr) data);
  return NULL;
}

/*****************************************************************************
*
*   SeqMgrClearFeatureIndexes clears every bioseq in an entity
*
*****************************************************************************/

static void SeqMgrClearIndexesProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  ObjMgrDataPtr  omdp;
  BoolPtr        rsult;

  if (sep == NULL || (! IS_Bioseq (sep))) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (DoSeqMgrFreeBioseqExtraData (omdp)) {
    rsult = (BoolPtr) mydata;
    *rsult = TRUE;
  }
}

NLM_EXTERN Boolean LIBCALL SeqMgrClearFeatureIndexes (Uint2 entityID, Pointer ptr)

{
  ObjMgrDataPtr  omdp;
  Boolean        rsult = FALSE;
  SeqEntryPtr    sep;

  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return FALSE;
  sep = SeqMgrGetTopSeqEntryForEntity (entityID);
  if (sep == NULL) return FALSE;
  SeqEntryExplore (sep, (Pointer) (&rsult), SeqMgrClearIndexesProc);

  /* clear out object manager time of indexing flag */

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL) {
    omdp->indexed = 0;
  }
  return rsult;
}

/*****************************************************************************
*
*   FindAppropriateBioseq finds the segmented bioseq if location is join on parts
*
*****************************************************************************/

static BioseqPtr FindAppropriateBioseq (SeqLocPtr loc, BioseqPtr tryfirst)

{
  BioseqPtr       bsp = NULL;
  BioseqExtraPtr  bspextra;
  BioseqSetPtr    bssp;
  ObjMgrDataPtr   omdp;
  BioseqPtr       part;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;
  SeqLocPtr       slp;

  if (loc == NULL) return NULL;
  sip = SeqLocId (loc);
  if (sip != NULL) {
    if (tryfirst != NULL && SeqIdIn (sip, tryfirst->id)) {
      bsp = tryfirst;
    } else {
      bsp = BioseqFind (sip);
    }

    /* first see if this is raw local part of segmented bioseq */

    if (bsp != NULL && bsp->repr == Seq_repr_raw) {
      omdp = SeqMgrGetOmdpForPointer (bsp);
      if (omdp != NULL && omdp->datatype == OBJ_BIOSEQ) {
        bspextra = (BioseqExtraPtr) omdp->extradata;
        if (bspextra != NULL) {
          if (bspextra->parentBioseq != NULL) {
            bsp = bspextra->parentBioseq;
          }
        }
      }
    }
    return bsp;
  }

  /* otherwise assume location is on multiple parts of a segmented set */

  slp = SeqLocFindNext (loc, NULL);
  if (slp == NULL) return NULL;
  sip = SeqLocId (slp);
  if (sip == NULL) return NULL;
  part = BioseqFind (sip);
  if (part == NULL) return NULL;
  omdp = SeqMgrGetOmdpForPointer (part);
  while (omdp != NULL) {
    if (omdp->datatype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) omdp->dataptr;
      if (bssp != NULL) {
        if (bssp->_class == BioseqseqSet_class_segset) {
          for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
            if (IS_Bioseq (sep)) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              if (bsp != NULL) {
                return bsp;
              }
            }
          }
        }
      }
    }
    omdp = SeqMgrGetOmdpForPointer (omdp->parentptr);
  }
  return NULL;
}

/*****************************************************************************
*
*   FindFirstLocalBioseq is called as a last resort if FindAppropriateBioseq
*     fails, and it scans the feature location to find the first local bioseq
*     referenced by a feature interval
*
*****************************************************************************/

static BioseqPtr FindFirstLocalBioseq (SeqLocPtr loc)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;
  SeqLocPtr  slp = NULL;

  if (loc == NULL) return NULL;

  while ((slp = SeqLocFindNext (loc, slp)) != NULL) {
    sip = SeqLocId (slp);
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL) return bsp;
    }
  }

  return NULL;
}

/*****************************************************************************
*
*   GetOffsetInFirstLocalBioseq is called to get the intervals on last resort bioseqs
*
*****************************************************************************/

static Int4 GetOffsetInFirstLocalBioseq (SeqLocPtr loc, BioseqPtr in, Uint1 which_end)

{
  SeqLocPtr  slp = NULL;
  Int4       val;

  if (loc == NULL) return -1;

  while ((slp = SeqLocFindNext (loc, slp)) != NULL) {
    val = GetOffsetInBioseq (slp, in, which_end);
    if (val != -1) return val;
  }

  return -1;
}

/*****************************************************************************
*
*   SeqMgrFindSMFeatItemPtr and SeqMgrFindSMFeatItemByID return SMFeatItemPtr
*     to access internal fields
*   SeqMgrGetDesiredDescriptor and SeqMgrGetDesiredFeature take an itemID,
*     position index, or SeqDescPtr or SeqFeatPtr, return the SeqDescPtr or
*     SeqFeatPtr, and fill in the context structure
*
*****************************************************************************/

NLM_EXTERN SMFeatItemPtr LIBCALL SeqMgrFindSMFeatItemPtr (SeqFeatPtr sfp)

{
  SMFeatItemPtr PNTR  array;
  BioseqPtr           bsp;
  BioseqExtraPtr      bspextra;
  SMFeatBlockPtr      curr;
  Int2                i;
  SMFeatItemPtr       item;
  Int4                L;
  Int4                mid;
  ObjMgrDataPtr       omdp;
  Int4                R;

  if (sfp == NULL) return NULL;
  bsp = FindAppropriateBioseq (sfp->location, NULL);
  if (bsp == NULL) {
    bsp = FindFirstLocalBioseq (sfp->location);
  }
  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;

  /* first try array sorted by SeqFeatPtr value */

  array = bspextra->featsBySfp;
  if (array != NULL && bspextra->numfeats > 0) {
    L = 0;
    R = bspextra->numfeats - 1;
    while (L < R) {
      mid = (L + R) / 2;
      item = array [mid];
      if (item != NULL && item->sfp < sfp) {
        L = mid + 1;
      } else {
        R = mid;
      }
    }

    item = array [R];
    if (item->sfp == sfp) return item;
  }

  /* now look in feature indices for cached feature information */

  curr = bspextra->featlisthead;
  while (curr != NULL) {

    if (curr->data != NULL) {
      for (i = 0; i < curr->index; i++) {
        item = &(curr->data [i]);
        if (item->sfp == sfp && (! item->ignore)) return item;
      }
    }

    curr = curr->next;
  }

  return NULL;
}

NLM_EXTERN SMFeatItemPtr LIBCALL SeqMgrFindSMFeatItemByID (BioseqPtr bsp, Uint2 itemID)

{
  SMFeatItemPtr PNTR  array;
  BioseqExtraPtr      bspextra;
  SMFeatBlockPtr      curr;
  Int2                i;
  SMFeatItemPtr       item;
  Int4                L;
  Int4                mid;
  ObjMgrDataPtr       omdp;
  Int4                R;

  if (bsp == NULL) return NULL;
  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;

  /* first try array sorted by itemID value */

  array = bspextra->featsByID;
  if (array != NULL && bspextra->numfeats > 0) {
    L = 0;
    R = bspextra->numfeats - 1;
    while (L < R) {
      mid = (L + R) / 2;
      item = array [mid];
      if (item != NULL && item->itemID < itemID) {
        L = mid + 1;
      } else {
        R = mid;
      }
    }

    item = array [R];
    if (item->itemID == itemID) return item;
  }

  /* now look in feature indices for cached feature information */

  curr = bspextra->featlisthead;
  while (curr != NULL) {

    if (curr->data != NULL) {
      for (i = 0; i < curr->index; i++) {
        item = &(curr->data [i]);
        if (item->itemID == itemID && (! item->ignore)) return item;
      }
    }

    curr = curr->next;
  }

  return NULL;
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetDesiredFeature (BioseqPtr bsp, Uint2 itemID,
                                                       Int2 index, SeqFeatPtr sfp,
                                                       SeqMgrFeatContext PNTR context)

{
  SMFeatItemPtr PNTR  array;
  BioseqExtraPtr      bspextra;
  SeqFeatPtr          curr;
  Uint2               entityID;
  SMFeatItemPtr       item = NULL;
  ObjMgrDataPtr       omdp;

  if (bsp == NULL) return NULL;
  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;

  if (itemID > 0) {
    item = SeqMgrFindSMFeatItemByID (bsp, itemID);
  } else if (index > 0) {
    array = bspextra->featsByPos;
    if (array != NULL && bspextra->numfeats > 0 && index <= bspextra->numfeats) {
      item = array [index - 1];
    }
  } else if (sfp != NULL) {
    item = SeqMgrFindSMFeatItemPtr (sfp);
  }
  if (item == NULL) return NULL;

  entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);

  curr = item->sfp;
  if (curr != NULL && context != NULL && (! item->ignore)) {
    context->entityID = entityID;
    context->itemID = item->itemID;
    context->sfp = curr;
    context->sap = item->sap;
    context->label = item->label;
    context->left = item->left;
    context->right = item->right;
    context->strand = item->strand;
    context->seqfeattype = FindFeatFromFeatDefType (item->subtype);;
    context->featdeftype = item->subtype;
    context->numivals = item->numivals;
    context->ivals = item->ivals;
    context->userdata = NULL;
    context->omdp = (Pointer) omdp;
    context->index = item->index + 1;
  }
  return curr;
}

NLM_EXTERN ValNodePtr LIBCALL SeqMgrGetDesiredDescriptor (BioseqPtr bsp, Uint2 itemID,
                                                          Int2 index, ValNodePtr sdp,
                                                          SeqMgrDescContext PNTR context)

{
  ValNodePtr         curr = NULL;
  SeqMgrDescContext  dfaultcontext;

  if (bsp == NULL) return NULL;

  if (context == NULL) {
    context = &dfaultcontext;
  }

  while ((curr = SeqMgrGetNextDescriptor (bsp, curr, 0, context)) != NULL) {
    if (itemID > 0 && itemID == context->itemID) return curr;
    if (index > 0 && index == context->index) return curr;
    if (sdp != NULL && sdp == curr) return curr;
  }

  return NULL;
}

/*****************************************************************************
*
*   RecordFeaturesInBioseqs callback explores bioseqs, bioseq sets, and features,
*     keeping a running total of the descriptor item counts, and records specific
*     information about features on each bioseq
*
*****************************************************************************/

typedef struct extraindex {
  BioseqPtr       lastbsp;
  SeqAnnotPtr     lastsap;
  BioseqSetPtr    lastbssp;
  SMSeqIdxPtr     segpartail;
  Int4            cumulative;
  Uint2           descrcount;
  Boolean         getLabel;
  Uint1           labeltype;
} ExtraIndex, PNTR ExtraIndexPtr;

static void SetDescriptorCounts (ValNodePtr sdp, ExtraIndexPtr exindx, Pointer thisitem)

{
  Uint2          count = 0;
  ObjMgrDataPtr  omdp;

  /* count bioseq or bioseq set descriptors, to calculate omdp.lastDescrItemID */

  if (sdp == NULL || exindx == NULL) return;
  omdp = SeqMgrGetOmdpForPointer (thisitem);
  if (omdp == NULL) return;

  omdp->lastDescrItemID = exindx->descrcount;
  while (sdp != NULL) {
    count++;
    sdp = sdp->next;
  }
  exindx->descrcount += count;
}

static void CreateBioseqExtraBlock (ObjMgrDataPtr omdp, BioseqPtr bsp)

{
  BioseqExtraPtr  bspextra;

  if (omdp == NULL || omdp->extradata != NULL) return;

  bspextra = (BioseqExtraPtr) MemNew (sizeof (BioseqExtra));
  omdp->extradata = (Pointer) bspextra;
  if (bspextra == NULL) return;

  omdp->reapextra = SeqMgrReapBioseqExtraFunc;
  omdp->reloadextra = SeqMgrReloadBioseqExtraFunc;
  omdp->freeextra = SeqMgrFreeBioseqExtraFunc;

  bspextra->bsp = bsp;
  bspextra->omdp = omdp;
  bspextra->min = INT4_MAX;
}

static SeqIdPtr SeqIdWithinBioseq (BioseqPtr bsp, SeqLocPtr slp)

{
  SeqIdPtr  a;
  SeqIdPtr  b;

  if (bsp == NULL || slp == NULL) return NULL;
  a = SeqLocId (slp);
  if (a == NULL) return NULL;
  for (b = bsp->id; b != NULL; b = b->next) {
    if (SeqIdComp (a, b) == SIC_YES) return b;
  }
  return NULL;
}

static void ProcessFeatureProducts (SeqFeatPtr sfp, Uint2 itemID)

{
  BioseqPtr       bsp;
  BioseqExtraPtr  bspextra;
  Int4            diff;
  Int4            min;
  ObjMgrDataPtr   omdp;
  SeqFeatPtr      prt;
  SeqAnnotPtr     sap;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  ValNode         vn;

  if (sfp == NULL || sfp->product == NULL) return;
  if (sfp->data.choice != SEQFEAT_CDREGION && sfp->data.choice != SEQFEAT_RNA) return;

  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) {
    CreateBioseqExtraBlock (omdp, bsp);
    bspextra = (BioseqExtraPtr) omdp->extradata;
  }
  if (bspextra == NULL) return;

  /* cds or rna reference stored in product bioseq's omdp.cdsOrRnaFeat */

  if (bspextra->cdsOrRnaFeat != NULL && bspextra->cdsOrRnaFeat != sfp) {
    ErrPostEx (SEV_WARNING, 0, 0, "SeqMgr indexing cds or rna progenitor already set");
  }
  if (omdp->tempload == TL_NOT_TEMP) {
    bspextra->cdsOrRnaFeat = sfp;
  }
  if (sfp->data.choice == SEQFEAT_RNA) return;

  /* if protFeat exists it was set by exhaustive gather on protein bioseq */

  if (bspextra->protFeat != NULL) return;

  /* calculate largest protein feature on cds's product bioseq */

  min = INT4_MAX;
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = (Pointer) bsp->id;
  vn.next = NULL;
  slp = (Pointer) (&vn);

  sap = bsp->annot;
  while (sap != NULL) {
    if (sap->type == 1) {
      prt = (SeqFeatPtr) sap->data;
      while (prt != NULL) {
        if (prt->data.choice == SEQFEAT_PROT) {

          /* get SeqId in bioseq that matches SeqId used for location */

          vn.data.ptrvalue = SeqIdWithinBioseq (bsp, prt->location);

          diff = SeqLocAinB (prt->location, slp);
          if (diff >= 0) {
            if (diff < min) {
              min = diff;
              if (omdp->tempload == TL_NOT_TEMP) {
                bspextra->protFeat = prt;
              }
            }
          }
        }
        prt = prt->next;
      }
    }
    sap = sap->next;
  }
}

static void RecordOneFeature (BioseqExtraPtr bspextra, ObjMgrDataPtr omdp,
                              BioseqPtr bsp, ExtraIndexPtr exindx, SeqFeatPtr sfp,
                              Int4 left, Int4 right, Uint2 itemID, Boolean ignore)

{
  Char            buf [129];
  SMFeatBlockPtr  curr;
  Int2            i;
  SMFeatItemPtr   item;
  Int4Ptr         ivals;
  SeqLocPtr       loc;
  SMFeatBlockPtr  next;
  Int2            numivals = 0;
  Boolean         single_interval;
  SeqLocPtr       slp = NULL;

  if (bspextra == NULL || omdp == NULL || bsp == NULL || exindx == NULL || sfp == NULL) return;

  if (bspextra->featlisttail != NULL) {

    /* just in case blocksize should was not set for some reason */

    if (bspextra->blocksize < 1) {
      bspextra->blocksize = 5;
    }

    curr = bspextra->featlisttail;
    if (curr->index >= bspextra->blocksize) {

      /* allocate next chunk in linked list of blocks */

      next = (SMFeatBlockPtr) MemNew (sizeof (SMFeatBlock));
      curr->next = next;

      if (next != NULL) {
        bspextra->featlisttail = next;
        curr = next;
      }
    }

    if (curr->index < bspextra->blocksize) {

      /* allocate data block if not yet done for this chunk */

      if (curr->data == NULL) {
        curr->data = (SMFeatItemPtr) MemNew (sizeof (SMFeatItem) * (size_t) (bspextra->blocksize));
      }

      /* now record desired information about current feature */

      if (curr->data != NULL) {
        item = &(curr->data [curr->index]);
        if (omdp->tempload == TL_NOT_TEMP) {
          item->sfp = sfp;
          item->sap = exindx->lastsap;
        }
        if (exindx->getLabel) {
          FeatDefLabel (sfp, buf, sizeof (buf) - 1, exindx->labeltype);
          item->label = StringSaveNoNull (buf);
        }
        item->left = left;
        item->right = right;
        item->strand = SeqLocStrand (sfp->location);
        item->subtype = FindFeatDefType (sfp);
        item->itemID = itemID;
        item->ignore = ignore;

        /* record start/stop pairs of intervals on target bioseq */

        single_interval = (Boolean) (item->subtype == FEATDEF_GENE ||
                                     item->subtype == FEATDEF_PUB);
        loc = SeqLocMerge (bsp, sfp->location, NULL, single_interval, TRUE, FALSE);

        slp = NULL;
        while ((slp = SeqLocFindNext (loc, slp)) != NULL) {
          numivals++;
        }
        if (numivals > 0) {
          ivals = MemNew (sizeof (Int4) * (numivals * 2));
          item->ivals = ivals;
          item->numivals = numivals;
          slp = NULL;
          i = 0;
          while ((slp = SeqLocFindNext (loc, slp)) != NULL) {
            ivals [i] = SeqLocStart (slp);
            i++;
            ivals [i] = SeqLocStop (slp);
            i++;
          }
        }
        SeqLocFree (loc);
      }

      /* increment count on current block */

      (curr->index)++;

      /* count all features */

      (bspextra->numfeats)++;

      /* count all gene, publication, and biosource features separately */

      if (sfp->data.choice == SEQFEAT_GENE) {
        (bspextra->numgenes)++;
      }
      if (sfp->data.choice == SEQFEAT_PUB) {
        (bspextra->numpubs)++;
      }
      if (sfp->data.choice == SEQFEAT_BIOSRC) {
        (bspextra->numorgs)++;
      }

    }
  }
}

/* callback for recording features and descriptor, prot, and cdsOrRna information */

static Boolean RecordFeaturesInBioseqs (GatherContextPtr gcp)

{
  BioseqPtr       bsp = NULL;
  BioseqExtraPtr  bspextra;
  BioseqSetPtr    bssp = NULL;
  Int2            count;
  Int4            diff;
  ExtraIndexPtr   exindx;
  Int4            left;
  ObjMgrDataPtr   omdp;
  Int4            right;
  SeqAnnotPtr     sap = NULL;
  ValNodePtr      sdp = NULL;
  SeqFeatPtr      sfp = NULL;
  SeqLocPtr       slp;
  SeqFeatPtr      tmp;
  Boolean         usingLocalBsp = FALSE;
  ValNode         vn;

  switch (gcp->thistype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) gcp->thisitem;
      if (bsp == NULL) return TRUE;
      sdp = bsp->descr;
      break;
    case OBJ_BIOSEQSET :
      bssp = (BioseqSetPtr) gcp->thisitem;
      if (bssp == NULL) return TRUE;
      sdp = bssp->descr;
      break;
    case OBJ_SEQANNOT :
      sap = (SeqAnnotPtr) gcp->thisitem;
      break;
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      break;
    default :
      return TRUE;
  }

  exindx = (ExtraIndexPtr) gcp->userdata;
  if (exindx == NULL) return FALSE;

  /* save bspItemID to support bioseq explore functions */

  if (bsp != NULL) {

    /* save last BioseqPtr to check first for appropriate bioseq */

    exindx->lastbsp = bsp;

    /* blocksize for new block based only on features packaged on bioseq */

    exindx->lastbssp = NULL;

    omdp = SeqMgrGetOmdpForPointer (bsp);
    if (omdp != NULL) {
      bspextra = (BioseqExtraPtr) omdp->extradata;
      if (bspextra == NULL) {
        CreateBioseqExtraBlock (omdp, bsp);
        bspextra = (BioseqExtraPtr) omdp->extradata;
      }
      if (bspextra != NULL) {
        bspextra->bspItemID = gcp->itemID;
      }
    }
  }

  /* save last BioseqSetPtr to calculate blocksize from bioseq set and bioseq features,
     features on bioseq set presumably being CDS or mRNA and applying only to nucleotides */

  if (bssp != NULL) {
    exindx->lastbssp = bssp;
  }

  /* count bioseq or bioseq set descriptors, to calculate lastDescrItemID */

  if (sdp != NULL) {
    SetDescriptorCounts (sdp, exindx, gcp->thisitem);
    return TRUE;
  }

  /* save SeqAnnotPtr containing next features to be gathered */

  if (sap != NULL) {
    exindx->lastsap = sap;
    return TRUE;
  }

  /* otherwise index features on every bioseq in entity */

  if (sfp == NULL) return TRUE;

  /* cds or rna reference stored in product bioseq's omdp.cdsOrRnaFeat,
     best protein feature in omdp.protFeat (do before adding CDS) */

  if (sfp->product != NULL) {
    ProcessFeatureProducts (sfp, gcp->itemID);
  }

  bsp = FindAppropriateBioseq (sfp->location, exindx->lastbsp);

  /* failure here can be due to SeqLoc that references far accession */

  if (bsp == NULL) {
    ErrPostItem (SEV_WARNING, 0, 0, "SeqMgr indexing feature location problem");

    /* if far accession, find first local bioseq on any location interval */

    bsp = FindFirstLocalBioseq (sfp->location);
    if (bsp == NULL) return TRUE;
    usingLocalBsp = TRUE;
  }

  /* assume subsequent features will be on this bioseq */

  exindx->lastbsp = bsp;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL) return TRUE;

  /* now prepare for adding feature to index */

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) {
    CreateBioseqExtraBlock (omdp, bsp);
    bspextra = (BioseqExtraPtr) omdp->extradata;
  }
  if (bspextra == NULL) return TRUE;

  /* get extreme left and right extents of feature location */

  if (usingLocalBsp) {

    left = GetOffsetInFirstLocalBioseq (sfp->location, bsp, SEQLOC_LEFT_END);
    if (left == -1) return TRUE;
    right = GetOffsetInFirstLocalBioseq (sfp->location, bsp, SEQLOC_RIGHT_END);
    if (right == -1) return TRUE;

  } else {

    left = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_LEFT_END);
    if (left == -1) return TRUE;
    right = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_RIGHT_END);
    if (right == -1) return TRUE;

  }

  /* if indexing protein bioseq, store largest protein feature */

  if (sfp->data.choice == SEQFEAT_PROT) {
    vn.choice = SEQLOC_WHOLE;
    vn.data.ptrvalue = (Pointer) bsp->id;
    vn.next = NULL;
    slp = (Pointer) &vn;

    /* get SeqId in bioseq that matches SeqId used for location */

    vn.data.ptrvalue = (Pointer) SeqIdWithinBioseq (bsp, sfp->location);

    diff = SeqLocAinB (sfp->location, slp);
    if (diff >= 0) {
      if (diff < bspextra->min) {
        bspextra->min = diff;
        if (omdp->tempload == TL_NOT_TEMP) {
          bspextra->protFeat = sfp;
        }
      }
    }
  }

  /* add feature item to linked list of blocks */

  if (bspextra->featlisthead == NULL) {
    bspextra->featlisthead = (SMFeatBlockPtr) MemNew (sizeof (SMFeatBlock));

    /* for first feature indexed on this bioseq, quickly see if few or many
       additional features, since most features on a bioseq are packaged in
       the same list, and most proteins only have one bioseq */

    for (tmp = sfp, count = 0;
         tmp != NULL && count < 50;
         tmp = tmp->next, count++) continue;

    /* extend count if above features were packaged on a bioseq set (presumably CDS or mRNA) */

    if (exindx->lastbssp != NULL) {
      for (sap = bsp->annot; sap != NULL; sap = sap->next) {
        if (sap->type == 1) {

          for (tmp = (SeqFeatPtr) sap->data;
               tmp != NULL && count < 50;
               tmp = tmp->next, count++) continue;

        }
      }
    }

    bspextra->blocksize = count;
  }
  if (bspextra->featlisttail == NULL) {
    bspextra->featlisttail = bspextra->featlisthead;
  }

  if (bspextra->featlisttail != NULL) {

    /* if feature spans origin, record with left < 0 */

    if (left > right && bsp->topology == TOPOLOGY_CIRCULAR) {
      left -= bsp->length;
    }

    RecordOneFeature (bspextra, omdp, bsp, exindx, sfp, left, right, gcp->itemID, FALSE);

    /* record gene, publication, and biosource features twice if spanning the origin */

    if (left < 0 && bsp->topology == TOPOLOGY_CIRCULAR) {
      if (sfp->data.choice == SEQFEAT_GENE ||
          sfp->data.choice == SEQFEAT_PUB ||
          sfp->data.choice == SEQFEAT_BIOSRC) {

        RecordOneFeature (bspextra, omdp, bsp, exindx, sfp, left + bsp->length,
                          right + bsp->length, gcp->itemID, TRUE);

      }
    }
  }

  return TRUE;
}

/*****************************************************************************
*
*   RecordSegmentsInBioseqs callback explores bioseq segments
*
*****************************************************************************/

static Boolean RecordSegmentsInBioseqs (GatherContextPtr gcp)

{
  BioseqPtr       bsp = NULL;
  BioseqExtraPtr  bspextra;
  Char            buf [80];
  ExtraIndexPtr   exindx;
  Int4            from;
  ObjMgrDataPtr   omdp;
  SMSeqIdxPtr     segpartptr;
  SeqIdPtr        sid;
  SeqIntPtr       sipp;
  SeqLocPtr       slp = NULL;
  Uint1           strand;
  Int4            to;

  switch (gcp->thistype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) gcp->thisitem;
      if (bsp == NULL) return TRUE;
      break;
    case OBJ_BIOSEQ_SEG :
      slp = (SeqLocPtr) gcp->thisitem;
      if (slp == NULL) return TRUE;
      break;
    default :
      return TRUE;
  }

  exindx = (ExtraIndexPtr) gcp->userdata;
  if (exindx == NULL) return FALSE;

  if (bsp != NULL) {
    if (bsp->repr == Seq_repr_seg) {
      exindx->lastbsp = bsp;
    } else {
      exindx->lastbsp = NULL;
    }
    exindx->cumulative = 0;
    return TRUE;
  }

  if (slp == NULL) return TRUE;

  bsp = exindx->lastbsp;
  if (bsp == NULL) return TRUE;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL) return TRUE;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) {
    CreateBioseqExtraBlock (omdp, bsp);
    bspextra = (BioseqExtraPtr) omdp->extradata;
  }
  if (bspextra == NULL) return TRUE;

  if (slp->choice == SEQLOC_INT && slp->data.ptrvalue != NULL) {
    sipp = (SeqIntPtr) (slp->data.ptrvalue);
    from = sipp->from;
    to = sipp->to;
    strand = sipp->strand;
  } else {
    from = 0;
    to = SeqLocLen (slp) - 1;
    strand = SeqLocStrand (slp);
  }

  if (to - from + 1 < 1) return TRUE;

  /* create and fill in SMSeqIdx element */

  segpartptr = MemNew (sizeof (SMSeqIdx));
  if (segpartptr != NULL) {
    sid = SeqLocId (slp);
    if (MakeReversedSeqIdString (sid, buf, sizeof (buf) - 1)) {
      segpartptr->slp = AsnIoMemCopy (slp,
                                      (AsnReadFunc) SeqLocAsnRead,
                                      (AsnWriteFunc) SeqLocAsnWrite);
      segpartptr->seqIdOfPart = StringSave (buf);
      segpartptr->parentBioseq = bsp;
      segpartptr->cumOffset = exindx->cumulative;
      segpartptr->from = from;
      segpartptr->to = to;
      segpartptr->strand = strand;
      segpartptr->itemID = gcp->itemID;
    }
  }

  exindx->cumulative += (to - from + 1);

  /* link into segparthead list of parts IDs */

  if (bspextra->segparthead == NULL) {
    bspextra->segparthead = segpartptr;
    exindx->segpartail = segpartptr;
  } else if (exindx->segpartail != NULL) {
    exindx->segpartail->next = segpartptr;
    exindx->segpartail = segpartptr;
  }

  return TRUE;
}

/*****************************************************************************
*
*   SortFeatItemListBySfp callback sorts array into feature item table by feature pointer
*   SortFeatItemListByPos callback sorts array into feature item table by feature position
*
*****************************************************************************/

static int LIBCALLBACK SortFeatItemListBySfp (VoidPtr vp1, VoidPtr vp2)

{
  SMFeatItemPtr PNTR  spp1 = vp1;
  SMFeatItemPtr PNTR  spp2 = vp2;
  SMFeatItemPtr       sp1;
  SMFeatItemPtr       sp2;

  if (spp1 == NULL || spp2 == NULL) return 0;
  sp1 = *((SMFeatItemPtr PNTR) spp1);
  sp2 = *((SMFeatItemPtr PNTR) spp2);
  if (sp1 == NULL || sp2 == NULL) return 0;

  if (sp1->sfp > sp2->sfp) {
    return 1;
  } else if (sp1->sfp < sp2->sfp) {
    return -1;
  }

  return 0;
}

static int LIBCALLBACK SortFeatItemListByPos (VoidPtr vp1, VoidPtr vp2)

{
  SMFeatItemPtr PNTR  spp1 = vp1;
  SMFeatItemPtr PNTR  spp2 = vp2;
  SMFeatItemPtr       sp1;
  SMFeatItemPtr       sp2;

  if (spp1 == NULL || spp2 == NULL) return 0;
  sp1 = *((SMFeatItemPtr PNTR) spp1);
  sp2 = *((SMFeatItemPtr PNTR) spp2);
  if (sp1 == NULL || sp2 == NULL) return 0;

  /* feature with smallest left extreme is first */

  if (sp1->left > sp2->left) {
    return 1;
  } else if (sp1->left < sp2->left) {
    return -1;

  /* reversing order so that longest feature is first */

  } else if (sp1->right > sp2->right) {
    return -1; /* was 1 */
  } else if (sp1->right < sp2->right) {
    return 1; /* was -1 */

  /* given identical extremes, put gene features first */

  } else if (sp1->subtype == FEATDEF_GENE) {
    return -1;
  } else if (sp2->subtype == FEATDEF_GENE) {
    return 1;

  /* then rna features */

  } else if (FindFeatFromFeatDefType (sp1->subtype) == SEQFEAT_RNA) {
    return -1;
  } else if (FindFeatFromFeatDefType (sp2->subtype) == SEQFEAT_RNA) {
    return 1;

  /* then cds features */

  } else if (sp1->subtype == FEATDEF_CDS) {
    return -1;
  } else if (sp2->subtype == FEATDEF_CDS) {
    return 1;
  }

  /* later will want to compare internal intervals */

  return 0;
}

/*****************************************************************************
*
*   IndexSegmentedParts callback builds index to speed up mapping
*     of parts to segmented bioseqs
*
*****************************************************************************/

static int LIBCALLBACK SortSeqIdxArray (VoidPtr ptr1, VoidPtr ptr2)

{
  Int2              compare;
  SMSeqIdxPtr PNTR  partp1 = ptr1;
  SMSeqIdxPtr PNTR  partp2 = ptr2;
  SMSeqIdxPtr       part1, part2;

  if (partp1 == NULL || partp2 == NULL) return 0;
  part1 = *((SMSeqIdxPtr PNTR) partp1);
  part2 = *((SMSeqIdxPtr PNTR) partp2);
  if (part1 == NULL || part2 == NULL) return 0;
  compare = StringCmp (part1->seqIdOfPart, part2->seqIdOfPart);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }
  return 0;
}

static void IndexSegmentedParts (SeqEntryPtr sep, BioseqPtr PNTR lastsegbsp)

{
  BioseqPtr         bsp;
  BioseqExtraPtr    bspextra;
  BioseqSetPtr      bssp;
  Int2              i;
  Int2              numsegs = 0;
  Int4              cumulative = 0;
  ObjMgrDataPtr     omdp;
  SMSeqIdxPtr PNTR  partsByLoc;
  SMSeqIdxPtr PNTR  partsBySeqId;
  SMSeqIdxPtr       segpartptr;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      IndexSegmentedParts (sep, lastsegbsp);
    }
    if (bssp->_class == BioseqseqSet_class_segset && lastsegbsp != NULL) {
      *lastsegbsp = NULL;
    }
    return;
  }

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;

  /* check for raw part packaged with segmented bioseq */

  if (bsp->repr == Seq_repr_raw && lastsegbsp != NULL && *lastsegbsp != NULL) {
    omdp = SeqMgrGetOmdpForPointer (bsp);
    if (omdp == NULL) return;

    bspextra = (BioseqExtraPtr) omdp->extradata;
    if (bspextra == NULL) {
      CreateBioseqExtraBlock (omdp, bsp);
      bspextra = (BioseqExtraPtr) omdp->extradata;
    }
    if (bspextra == NULL) return;

    /* now record segmented parent of raw part if all are packaged together */

    bspextra->parentBioseq = *lastsegbsp;
    return;
  }

  if (bsp->repr != Seq_repr_seg) return;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL) return;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) {
    CreateBioseqExtraBlock (omdp, bsp);
    bspextra = (BioseqExtraPtr) omdp->extradata;
  }
  if (bspextra == NULL) return;

  if (lastsegbsp != NULL) {
    *lastsegbsp = bsp;
  }

  for (segpartptr = bspextra->segparthead;
       segpartptr != NULL;
       segpartptr = segpartptr->next) {
    numsegs++;
  }

  bspextra->numsegs = numsegs;
  segpartptr = bspextra->segparthead;
  if (numsegs < 1 || segpartptr == NULL) return;

  partsByLoc = (SMSeqIdxPtr PNTR) MemNew (sizeof (SMSeqIdxPtr) * (numsegs + 1));
  bspextra->partsByLoc = partsByLoc;

  if (partsByLoc != NULL) {
    i = 0;
    while (i < numsegs && segpartptr != NULL) {
      partsByLoc [i] = segpartptr;
      segpartptr = segpartptr->next;
      i++;
    }

    partsBySeqId = (SMSeqIdxPtr PNTR) MemNew (sizeof (SMSeqIdxPtr) * (numsegs + 1));
    bspextra->partsBySeqId = partsBySeqId;

    if (partsBySeqId != NULL) {
      for (i = 0; i < numsegs; i++) {
        partsBySeqId [i] = partsByLoc [i];
      }

      /* sort array by SeqId for binary search */

      HeapSort ((Pointer) partsBySeqId, numsegs, sizeof (SMSeqIdxPtr), SortSeqIdxArray);
    }

  }
}

/*****************************************************************************
*
*   IndexRecordedFeatures callback builds sorted arrays of features and genes
*
*****************************************************************************/

static void IndexRecordedFeatures (SeqEntryPtr sep)

{
  BioseqPtr           bsp;
  BioseqExtraPtr      bspextra;
  BioseqSetPtr        bssp;
  SMFeatBlockPtr      curr;
  SMFeatItemPtr PNTR  featsByID;
  SMFeatItemPtr PNTR  featsBySfp;
  SMFeatItemPtr PNTR  featsByPos;
  SMFeatItemPtr PNTR  genesByPos;
  SMFeatItemPtr PNTR  pubsByPos;
  SMFeatItemPtr PNTR  orgsByPos;
  Int4                i;
  Int4                j;
  SMFeatItemPtr       item;
  Int4                numfeats;
  Int4                numgenes;
  Int4                numpubs;
  Int4                numorgs;
  ObjMgrDataPtr       omdp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      IndexRecordedFeatures (sep);
    }
    return;
  }

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL) return;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return;

  numfeats = bspextra->numfeats;
  numgenes = bspextra->numgenes;
  numpubs = bspextra->numpubs;
  numorgs = bspextra->numorgs;

  curr = bspextra->featlisthead;

  if (bspextra->numfeats > 0 && curr != NULL) {

    /* build array of pointers into feature items */

    featsByID = (SMFeatItemPtr PNTR) MemNew (sizeof (SMFeatItemPtr) * (numfeats + 1));
    bspextra->featsByID = featsByID;

    if (featsByID != NULL) {
      i = 0;
      j = 0;
      while (i < numfeats && curr != NULL) {
        if (j >= curr->index || j >= bspextra->blocksize) {
          j = 0;
          curr = curr->next;
        }
        if (curr != NULL && j < curr->index && curr->data != NULL) {
          featsByID [i] = &(curr->data [j]);
          i++;
          j++;
        }
      }
      if (i < numfeats) {
        ErrPostEx (SEV_WARNING, 0, 0, "SeqMgr indexing feature table build problem");
      }

      featsBySfp = (SMFeatItemPtr PNTR) MemNew (sizeof (SMFeatItemPtr) * (numfeats + 1));
      bspextra->featsBySfp = featsBySfp;

      if (featsBySfp != NULL) {
        for (i = 0; i < numfeats; i++) {
          featsBySfp [i] = featsByID [i];
        }

        /* sort all features by SeqFeatPtr value */

        HeapSort ((VoidPtr) featsBySfp, (size_t) numfeats, sizeof (SMFeatItemPtr), SortFeatItemListBySfp);
      }

      featsByPos = (SMFeatItemPtr PNTR) MemNew (sizeof (SMFeatItemPtr) * (numfeats + 1));
      bspextra->featsByPos = featsByPos;

      if (featsByPos != NULL) {
        for (i = 0; i < numfeats; i++) {
          featsByPos [i] = featsByID [i];
        }

        /* sort all features by feature location on bioseq */

        HeapSort ((VoidPtr) featsByPos, (size_t) numfeats, sizeof (SMFeatItemPtr), SortFeatItemListByPos);

        for (i = 0; i < numfeats; i++) {
          item = featsByPos [i];
          if (item != NULL) {
            item->index = i;
          }
        }

        /* build subarray of sorted gene features for lookup by overlap */

        if (numgenes > 0) {

          genesByPos = (SMFeatItemPtr PNTR) MemNew (sizeof (SMFeatItemPtr) * (numgenes + 1));
          bspextra->genesByPos = genesByPos;

          if (genesByPos != NULL) {
            i = 0;
            j = 0;
            while (i < numfeats && j < numgenes) {
              item = featsByPos [i];
              if (item->subtype == FEATDEF_GENE) {
                genesByPos [j] = item;
                j++;
              }
              i++;
            }
          }
        }

        /* build subarray of sorted publication features for lookup by overlap */

        if (numpubs > 0) {

          pubsByPos = (SMFeatItemPtr PNTR) MemNew (sizeof (SMFeatItemPtr) * (numpubs + 1));
          bspextra->pubsByPos = pubsByPos;

          if (pubsByPos != NULL) {
            i = 0;
            j = 0;
            while (i < numfeats && j < numpubs) {
              item = featsByPos [i];
              if (item->subtype == FEATDEF_PUB) {
                pubsByPos [j] = item;
                j++;
              }
              i++;
            }
          }
        }

        /* build subarray of sorted biosource features for lookup by overlap */

        if (numorgs > 0) {

          orgsByPos = (SMFeatItemPtr PNTR) MemNew (sizeof (SMFeatItemPtr) * (numorgs + 1));
          bspextra->orgsByPos = orgsByPos;

          if (orgsByPos != NULL) {
            i = 0;
            j = 0;
            while (i < numfeats && j < numorgs) {
              item = featsByPos [i];
              if (item->subtype == FEATDEF_BIOSRC) {
                orgsByPos [j] = item;
                j++;
              }
              i++;
            }
          }
        }
      }

    }
  }
}

/*****************************************************************************
*
*   SeqMgrReindexBioseqExtraData refreshes internal indices for rapid retrieval
*
*****************************************************************************/

NLM_EXTERN Uint2 LIBCALL SeqMgrIndexFeatures (Uint2 entityID, Pointer ptr, Uint1 labeltype)

{
  ExtraIndex     exind;
  GatherScope    gs;
  BioseqPtr      lastsegbsp = NULL;
  SeqEntryPtr    oldscope;
  ObjMgrDataPtr  omdp;
  SeqEntryPtr    sep;

  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return 0;

  /* reset any existing index data on all bioseqs in entity */

  SeqMgrClearFeatureIndexes (entityID, NULL);

  /* want to scope to bioseqs within the entity, to allow for colliding IDs */

  sep = SeqMgrGetTopSeqEntryForEntity (entityID);

  /* set scope for FindAppropriateBioseq, FindFirstLocalBioseq */

  oldscope = SeqEntrySetScope (sep);

  /* gather all segmented locations */

  exind.lastbsp = NULL;
  exind.lastsap = NULL;
  exind.lastbssp = NULL;
  exind.segpartail = NULL;
  exind.descrcount = 0;
  exind.getLabel = TRUE;
  exind.labeltype = labeltype;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  MemSet ((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore [OBJ_BIOSEQ] = FALSE;
  gs.ignore [OBJ_BIOSEQSET] = FALSE;
  gs.ignore [OBJ_BIOSEQ_SEG] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) (&exind), RecordSegmentsInBioseqs, &gs);

  /* build indexes to speed mapping of parts to segmented bioseq */

  lastsegbsp = NULL;

  IndexSegmentedParts (sep, &lastsegbsp);

  /* now gather to get descriptor itemID counts on each bioseq or bioseq set,
     and record features on the bioseq indicated by the feature location */

  exind.lastbsp = NULL;
  exind.lastsap = NULL;
  exind.lastbssp = NULL;
  exind.segpartail = NULL;
  exind.descrcount = 0;
  exind.getLabel = TRUE;
  exind.labeltype = labeltype;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  MemSet ((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore [OBJ_BIOSEQ] = FALSE;
  gs.ignore [OBJ_BIOSEQSET] = FALSE;
  gs.ignore [OBJ_SEQANNOT] = FALSE;
  gs.ignore [OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) (&exind), RecordFeaturesInBioseqs, &gs);

  /* finish building array of sorted features on each indexed bioseq */

  IndexRecordedFeatures (sep);

  /* resetset scope used to limit FindAppropriateBioseq, FindFirstLocalBioseq */

  SeqEntrySetScope (oldscope);

  /* stamp top of entity with time of indexing */

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL) {
    omdp->indexed = GetSecs ();
  }

  return entityID;
}

/*****************************************************************************
*
*   SeqMgrIsBioseqIndexed checks for presence of time of indexing stamp
*
*****************************************************************************/

NLM_EXTERN time_t LIBCALL SeqMgrFeaturesAreIndexed (Uint2 entityID)

{
  ObjMgrDataPtr  omdp;

  if (entityID == 0) return 0;
  omdp = ObjMgrGetData (entityID);
  if (omdp == NULL) return 0;
  return omdp->indexed;
}

/*****************************************************************************
*
*   SeqMgrGetBestProteinFeature and SeqMgrGetCDSgivenProduct take a protein
*     bioseq to get the best protein feature or encoding CDS
*   SeqMgrGetRNAgivenProduct takes an mRNA (cDNA) bioseq and gets encoding mRNA
*     feature on the genomic bioseq
*
*****************************************************************************/

static void SetContextForFeature (SeqFeatPtr sfp, SeqMgrFeatContext PNTR context, ObjMgrDataPtr omdp)

{
  SMFeatItemPtr  best;

  if (sfp == NULL || context == NULL || omdp == NULL) return;
  best = SeqMgrFindSMFeatItemPtr (sfp);
  if (best == NULL) return;
  context->entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);
  context->itemID = best->itemID;
  context->sfp = best->sfp;
  context->sap = best->sap;
  context->label = best->label;
  context->left = best->left;
  context->right = best->right;
  context->strand = best->strand;
  context->seqfeattype = FindFeatFromFeatDefType (best->subtype);
  context->featdeftype = best->subtype;
  context->numivals = best->numivals;
  context->ivals = best->ivals;
  context->userdata = NULL;
  context->omdp = (Pointer) omdp;
  context->index = best->index + 1;
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetBestProteinFeature (BioseqPtr bsp,
                                                           SeqMgrFeatContext PNTR context)

{
  BioseqExtraPtr  bspextra;
  ObjMgrDataPtr   omdp;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;
  SetContextForFeature (bspextra->protFeat, context, omdp);
  return bspextra->protFeat;
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetCDSgivenProduct (BioseqPtr bsp,
                                                        SeqMgrFeatContext PNTR context)

{
  BioseqExtraPtr  bspextra;
  ObjMgrDataPtr   omdp;
  SeqFeatPtr      sfp;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;
  sfp = bspextra->cdsOrRnaFeat;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return NULL;
  SetContextForFeature (sfp, context, omdp);
  return sfp;
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetRNAgivenProduct (BioseqPtr bsp,
                                                        SeqMgrFeatContext PNTR context)

{
  BioseqExtraPtr  bspextra;
  ObjMgrDataPtr   omdp;
  SeqFeatPtr      sfp;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;
  sfp = bspextra->cdsOrRnaFeat;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return NULL;
  SetContextForFeature (sfp, context, omdp);
  return sfp;
}

/*****************************************************************************
*
*   SeqMgrGetGeneXref, SeqMgrGeneIsSuppressed, SeqMgrGetOverlappingGene,
*     and SeqMgrGetOverlappingPub
*
*****************************************************************************/

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

NLM_EXTERN GeneRefPtr LIBCALL SeqMgrGetGeneXref (SeqFeatPtr sfp)

{
  GeneRefPtr      grp = NULL;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return NULL;
  xref = sfp->xref;
  while (xref != NULL && xref->data.choice != SEQFEAT_GENE) {
    xref = xref->next;
  }
  if (xref != NULL) {
    grp = (GeneRefPtr) xref->data.value.ptrvalue;
  }
  return grp;
}

NLM_EXTERN Boolean LIBCALL SeqMgrGeneIsSuppressed (GeneRefPtr grp)

{
  if (grp == NULL) return FALSE;
  if (grp != NULL && HasNoText (grp->locus) && HasNoText (grp->allele) &&
      HasNoText (grp->desc) && HasNoText (grp->maploc) &&
      grp->db == NULL && grp->syn == NULL) return TRUE;
  return FALSE;
}

static SeqFeatPtr SeqMgrGetBestOverlappingFeat (SeqFeatPtr sfp, Uint2 subtype, SeqMgrFeatContext PNTR context)

{
  SMFeatItemPtr PNTR  array = NULL;
  SMFeatItemPtr       best = NULL;
  BioseqPtr           bsp;
  BioseqExtraPtr      bspextra;
  Int4                diff;
  SMFeatItemPtr       feat;
  Int2                index = 0;
  Int4                L;
  Int4                left;
  Int4                max;
  Int4                mid;
  Int4                num = 0;
  ObjMgrDataPtr       omdp;
  Int4                R;
  Int4                right;
  Uint1               strand;

  if (sfp == NULL) return NULL;
  bsp = FindAppropriateBioseq (sfp->location, NULL);
  if (bsp == NULL) {
    bsp = FindFirstLocalBioseq (sfp->location);
  }
  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;
  switch (subtype) {
    case FEATDEF_GENE :
      array = bspextra->genesByPos;
      num = bspextra->numgenes;
      break;
    case FEATDEF_PUB :
      array = bspextra->pubsByPos;
      num = bspextra->numpubs;
      break;
    case FEATDEF_BIOSRC :
      array = bspextra->orgsByPos;
      num = bspextra->numorgs;
      break;
    default :
      break;
  }
  if (array == NULL || num < 1) return NULL;

  left = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_LEFT_END);
  if (left == -1) return NULL;
  right = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_RIGHT_END);
  if (right == -1) return NULL;

  /* if feature spans origin, normalize with left < 0 */

  if (left > right && bsp->topology == TOPOLOGY_CIRCULAR) {
    left -= bsp->length;
  }

  /* binary search to leftmost candidate within the genesByPos, pubsByPos, or orgsByPos array */

  L = 0;
  R = num - 1;
  while (L < R) {
    mid = (L + R) / 2;
    feat = array [mid];
    if (feat != NULL && feat->right < left) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  /* linear scan to smallest covering gene, publication, or biosource */

  best = NULL;
  index = 0;

  feat = array [R];
  max = INT4_MAX;
  strand = SeqLocStrand (sfp->location);
  while (R < num && feat != NULL && feat->left <= right) {
    if (feat->left <= left && feat->right >= right) {
      if (feat->strand == strand ||
          strand == Seq_strand_unknown ||
          feat->strand == Seq_strand_unknown) {
        diff = (left - feat->left) + (feat->right - right);
        if (diff < max) {
          best = feat;
          index = R;
          max = diff;
        }
      }
    }
    R++;
    feat = array [R];
  }

  if (best != NULL) {
    if (context != NULL) {
      context->entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);
      context->itemID = best->itemID;
      context->sfp = best->sfp;
      context->sap = best->sap;
      context->label = best->label;
      context->left = best->left;
      context->right = best->right;
      context->strand = best->strand;
      context->seqfeattype = FindFeatFromFeatDefType (best->subtype);
      context->featdeftype = best->subtype;
      context->numivals = best->numivals;
      context->ivals = best->ivals;
      context->userdata = NULL;
      context->omdp = (Pointer) omdp;
      context->index = best->index + 1;
    }
    return best->sfp;
  }

  return NULL;
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetOverlappingGene (SeqFeatPtr sfp, SeqMgrFeatContext PNTR context)

{
  return SeqMgrGetBestOverlappingFeat (sfp, FEATDEF_GENE, context);
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetOverlappingPub (SeqFeatPtr sfp, SeqMgrFeatContext PNTR context)

{
  return SeqMgrGetBestOverlappingFeat (sfp, FEATDEF_PUB, context);
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetOverlappingSource (SeqFeatPtr sfp, SeqMgrFeatContext PNTR context)

{
  return SeqMgrGetBestOverlappingFeat (sfp, FEATDEF_BIOSRC, context);
}

/*****************************************************************************
*
*   SeqMgrGetNextDescriptor and SeqMgrGetNextFeature
*
*****************************************************************************/

NLM_EXTERN ValNodePtr LIBCALL SeqMgrGetNextDescriptor (BioseqPtr bsp, ValNodePtr curr,
                                                       Uint1 seqDescChoice,
                                                       SeqMgrDescContext PNTR context)

{
  BioseqSetPtr   bssp;
  Uint2          entityID;
  ObjMgrDataPtr  omdp;
  SeqEntryPtr    sep;
  ValNode        vn;

  if (context == NULL) return NULL;

  /* if curr is NULL, initialize context fields (in user's stack) */

  if (curr == NULL) {
    if (bsp == NULL) return NULL;
    omdp = SeqMgrGetOmdpForPointer (bsp);
    if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;

    context->omdp = (Pointer) omdp;
    context->itemID = omdp->lastDescrItemID;
    context->index = 0;

    /* start curr just before beginning of bioseq descriptor list */

    curr = &vn;
    vn.choice = 0;
    vn.data.ptrvalue = 0;
    vn.next = bsp->descr;
  }

  omdp = (ObjMgrDataPtr) context->omdp;
  if (omdp == NULL) return NULL;

  entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);
  sep = ObjMgrGetChoiceForData (omdp->dataptr);

  /* now look for next appropriate descriptor after curr in current chain */

  while (curr != NULL) {
    curr = curr->next;
    if (curr != NULL) {
      (context->itemID)++;
      (context->index)++;
      if (seqDescChoice == 0 || curr->choice == seqDescChoice) {
        context->entityID = entityID;
        context->sdp = curr;
        context->sep = sep;
        context->seqdesctype = curr->choice;
        context->userdata = NULL;
        context->omdp = (Pointer) omdp;
        return curr;
      }
    }
  }

  /* now go up omdp chain looking for next descriptor */

  while (curr == NULL) {
    omdp = SeqMgrGetOmdpForPointer (omdp->parentptr);
    if (omdp == NULL) return NULL;

    /* update current omdp in context */

    context->omdp = (Pointer) omdp;
    context->itemID = omdp->lastDescrItemID;

    switch (omdp->datatype) {
      case OBJ_BIOSEQ :
        bsp = (BioseqPtr) omdp->dataptr;
        curr = bsp->descr;
        break;
      case OBJ_BIOSEQSET :
        bssp = (BioseqSetPtr) omdp->dataptr;
        curr = bssp->descr;
        break;
      default :
        break;
    }

    sep = ObjMgrGetChoiceForData (omdp->dataptr);

    /* now look for first appropriate descriptor in current chain */

    while (curr != NULL) {
      (context->itemID)++;
      (context->index)++;
      if (seqDescChoice == 0 || curr->choice == seqDescChoice) {
        context->entityID = entityID;
        context->sdp = curr;
        context->sep = sep;
        context->seqdesctype = curr->choice;
        context->userdata = NULL;
        context->omdp = (Pointer) omdp;
        return curr;
      }
      curr = curr->next;
    }
  }

  return curr;
}

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetNextFeature (BioseqPtr bsp, SeqFeatPtr curr,
                                                    Uint1 seqFeatChoice, Uint1 featDefChoice,
                                                    SeqMgrFeatContext PNTR context)

{
  BioseqExtraPtr      bspextra;
  Uint2               entityID;
  SMFeatItemPtr PNTR  featsByPos;
  Int2                i;
  SMFeatItemPtr       item;
  ObjMgrDataPtr       omdp;
  Uint1               seqfeattype;

  if (context == NULL) return NULL;

  /* if curr is NULL, initialize context fields (in user's stack) */

  if (curr == NULL) {
    if (bsp == NULL) return NULL;
    omdp = SeqMgrGetOmdpForPointer (bsp);
    if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;

    context->omdp = (Pointer) omdp;
    context->index = 0;
  }

  omdp = (ObjMgrDataPtr) context->omdp;
  if (omdp == NULL) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;
  featsByPos = bspextra->featsByPos;
  if (featsByPos == NULL || bspextra->numfeats < 1) return NULL;

  entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);

  i = context->index;

  /* now look for next appropriate feature */

  while (i < bspextra->numfeats) {
    item = featsByPos [i];
    if (item != NULL) {
      curr = item->sfp;
      i++;
      if (curr != NULL) {
        seqfeattype = FindFeatFromFeatDefType (item->subtype);
        if ((seqFeatChoice == 0 || seqfeattype == seqFeatChoice) &&
            (featDefChoice == 0 || item->subtype == featDefChoice) &&
            (! item->ignore)) {
          context->entityID = entityID;
          context->itemID = item->itemID;
          context->sfp = curr;
          context->sap = item->sap;
          context->label = item->label;
          context->left = item->left;
          context->right = item->right;
          context->strand = item->strand;
          context->seqfeattype = seqfeattype;
          context->featdeftype = item->subtype;
          context->numivals = item->numivals;
          context->ivals = item->ivals;
          context->userdata = NULL;
          context->omdp = (Pointer) omdp;
          context->index = item->index + 1;
          return curr;
        }
      }
    }
  }

  return NULL;
}

/*****************************************************************************
*
*   SeqMgrExploreBioseqs, SeqMgrExploreDescriptors, and SeqMgrExploreFeatures
*
*****************************************************************************/

static Boolean JustExamineBioseqs (SeqEntryPtr sep, BioseqSetPtr bssp,
                                   SeqMgrBioseqContextPtr context,
                                   SeqMgrBioseqExploreProc userfunc,
                                   Boolean nucs, Boolean prots, Boolean parts)

{
  BioseqPtr       bsp;
  BioseqExtraPtr  bspextra;
  ObjMgrDataPtr   omdp;

  if (sep == NULL || context == NULL || userfunc == NULL) return FALSE;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return TRUE;

    /* check for desired molecule type */

    if (ISA_na (bsp->mol) && (! nucs)) return TRUE;
    if (ISA_aa (bsp->mol) && (! prots)) return TRUE;

    omdp = SeqMgrGetOmdpForPointer (bsp);
    if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return TRUE;
    bspextra = (BioseqExtraPtr) omdp->extradata;
    if (bspextra == NULL) return TRUE;

    context->itemID = bspextra->bspItemID;
    context->bsp = bsp;
    context->sep = sep;
    context->bssp = bssp;
    context->omdp = omdp;
    (context->index)++;

    /* continue until user function returns FALSE, then exit all recursions */

    if (! userfunc (bsp, context)) return FALSE;
    return TRUE;
  }

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return TRUE;

    /* check to see if parts should be explored */

    if (bssp->_class == BioseqseqSet_class_parts && (! parts)) return TRUE;

    /* recursively explore bioseq set until user function returns FALSE */

    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      if (! JustExamineBioseqs (sep, bssp, context, userfunc, nucs, prots, parts)) return FALSE;
    }
  }

  return TRUE;
}

NLM_EXTERN Boolean LIBCALL SeqMgrExploreBioseqs (Uint2 entityID, Pointer ptr, Pointer userdata,
                                                 SeqMgrBioseqExploreProc userfunc,
                                                 Boolean nucs, Boolean prots, Boolean parts)

{
  SeqMgrBioseqContext  context;
  SeqEntryPtr          sep;

  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return 0;
  sep = SeqMgrGetTopSeqEntryForEntity (entityID);
  if (sep == NULL) return FALSE;
  if (userfunc == NULL) return FALSE;

  context.entityID = entityID;
  context.index = 0;
  context.userdata = userdata;

  /* recursive call to explore SeqEntry and pass appropriate bioseqs to user */

  JustExamineBioseqs (sep, NULL, &context, userfunc, nucs, prots, parts);

  return TRUE;
}

NLM_EXTERN Boolean LIBCALL SeqMgrExploreSegments (BioseqPtr bsp, Pointer userdata,
                                                  SeqMgrSegmentExploreProc userfunc)

{
  BioseqExtraPtr        bspextra;
  SeqMgrSegmentContext  context;
  Uint2                 entityID;
  Int2                  i;
  ObjMgrDataPtr         omdp;
  SMSeqIdxPtr PNTR      partsByLoc;
  SMSeqIdxPtr           segpartptr;
  SeqLocPtr             slp;

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return FALSE;
  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return FALSE;
  if (userfunc == NULL) return FALSE;
  entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return FALSE;
  partsByLoc = bspextra->partsByLoc;
  if (partsByLoc == NULL || bspextra->numsegs < 1) return FALSE;

  for (i = 0; i < bspextra->numsegs; i++) {
    segpartptr = partsByLoc [i];
    if (segpartptr != NULL) {
      slp = segpartptr->slp;
      context.entityID = entityID;
      context.itemID = segpartptr->itemID;
      context.slp = slp;
      context.parent = segpartptr->parentBioseq;
      context.cumOffset = segpartptr->cumOffset;
      context.from = segpartptr->from;
      context.to = segpartptr->to;
      context.strand = segpartptr->strand;
      context.userdata = userdata;
      context.omdp = (Pointer) omdp;
      context.index = i + 1;
      if (! userfunc (slp, &context)) return TRUE;
    }
  }

  return TRUE;
}

NLM_EXTERN Boolean LIBCALL SeqMgrExploreDescriptors (BioseqPtr bsp, Pointer userdata,
                                                     SeqMgrDescExploreProc userfunc,
                                                     BoolPtr seqDescFilter)

{
  BioseqSetPtr       bssp;
  SeqMgrDescContext  context;
  Uint2              entityID;
  Uint2              itemID;
  ObjMgrDataPtr      omdp;
  ValNodePtr         sdp;
  SeqEntryPtr        sep;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return FALSE;
  if (userfunc == NULL) return FALSE;
  entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);

  context.index = 0;
  while (omdp != NULL) {
    itemID = omdp->lastDescrItemID;
    sdp = NULL;
    switch (omdp->datatype) {
      case OBJ_BIOSEQ :
        bsp = (BioseqPtr) omdp->dataptr;
        sdp = bsp->descr;
        break;
      case OBJ_BIOSEQSET :
        bssp = (BioseqSetPtr) omdp->dataptr;
        sdp = bssp->descr;
        break;
      default :
        break;
    }

    sep = ObjMgrGetChoiceForData (omdp->dataptr);

    /* call for every appropriate descriptor in current chain */

    while (sdp != NULL) {
      itemID++;
      if (seqDescFilter == NULL || seqDescFilter [sdp->choice]) {
        context.entityID = entityID;
        context.itemID = itemID;
        context.sdp = sdp;
        context.sep = sep;
        context.seqdesctype = sdp->choice;
        context.userdata = userdata;
        context.omdp = (Pointer) omdp;
        (context.index)++;
        if (! userfunc (sdp, &context)) return TRUE;
      }
      sdp = sdp->next;
    }

    /* now go up omdp chain looking for next descriptor */

    omdp = SeqMgrGetOmdpForPointer (omdp->parentptr);
  }
  return TRUE;
}

NLM_EXTERN Boolean LIBCALL SeqMgrExploreFeatures (BioseqPtr bsp, Pointer userdata,
                                                  SeqMgrFeatExploreProc userfunc,
                                                  SeqLocPtr locationFilter,
                                                  BoolPtr seqFeatFilter, BoolPtr featDefFilter)

{
  BioseqExtraPtr      bspextra;
  SeqMgrFeatContext   context;
  Uint2               entityID;
  SMFeatItemPtr PNTR  featsByPos;
  Int2                i;
  SMFeatItemPtr       item;
  Int4                left = INT4_MIN;
  ObjMgrDataPtr       omdp;
  Int4                right = INT4_MAX;
  Uint1               seqfeattype;
  SeqFeatPtr          sfp;

  omdp = SeqMgrGetOmdpForPointer (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return FALSE;
  if (userfunc == NULL) return FALSE;
  entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return FALSE;
  featsByPos = bspextra->featsByPos;
  if (featsByPos == NULL || bspextra->numfeats < 1) return FALSE;

  if (locationFilter != NULL) {
    left = GetOffsetInBioseq (locationFilter, bsp, SEQLOC_LEFT_END);
    if (left == -1) left = INT4_MIN;
    right = GetOffsetInBioseq (locationFilter, bsp, SEQLOC_RIGHT_END);
    if (right == -1) right = INT4_MAX;
  }

  /* call for every appropriate feature in sorted list */

  for (i = 0; i < bspextra->numfeats; i++) {
    item = featsByPos [i];
    if (item != NULL) {
      sfp = item->sfp;
      seqfeattype = FindFeatFromFeatDefType (item->subtype);
      if ((seqFeatFilter == NULL || seqFeatFilter [seqfeattype]) &&
          (featDefFilter == NULL || featDefFilter [item->subtype]) &&
          (locationFilter == NULL || (item->right >= left && item->left <= right)) &&
          (! item->ignore)) {
        context.entityID = entityID;
        context.itemID = item->itemID;
        context.sfp = sfp;
        context.sap = item->sap;
        context.label = item->label;
        context.left = item->left;
        context.right = item->right;
        context.strand = item->strand;
        context.seqfeattype = seqfeattype;
        context.featdeftype = item->subtype;
        context.numivals = item->numivals;
        context.ivals = item->ivals;
        context.userdata = userdata;
        context.omdp = (Pointer) omdp;
        context.index = item->index + 1;
        if (! userfunc (sfp, &context)) return TRUE;
      }
    }
  }
  return TRUE;
}

/*****************************************************************************
*
*   SeqMgrMapPartToSegmentedBioseq can speed up sequtil's CheckPointInBioseq
*     for indexed part bioseq to segmented bioseq mapping
*
*****************************************************************************/

static SMSeqIdxPtr BinarySearchPartToSegmentMap (BioseqPtr in, Int4 pos, BioseqPtr bsp, SeqIdPtr sip)

{
  BioseqExtraPtr    bspextra;
  Char              buf [80];
  Int2              compare;
  ObjMgrDataPtr     omdp;
  SMSeqIdxPtr PNTR  partsBySeqId;
  SMSeqIdxPtr       segpartptr;
  Int2              L, R, mid;

  if (in == NULL) return NULL;
  omdp = SeqMgrGetOmdpForPointer (in);
  if (omdp == NULL) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;

  partsBySeqId = bspextra->partsBySeqId;
  if (partsBySeqId == NULL || bspextra->numsegs < 1) return NULL;

  if (bsp != NULL) {
    sip = bsp->id;
  }
  if (sip == NULL) return NULL;

  /* binary search into array on segmented bioseq sorted by part seqID (reversed) string */

  while (sip != NULL) {
    if (MakeReversedSeqIdString (sip, buf, sizeof (buf) - 1)) {
      L = 0;
      R = bspextra->numsegs - 1;
      while (L < R) {
        mid = (L + R) / 2;
        segpartptr = partsBySeqId [mid];
        compare = StringCmp (segpartptr->seqIdOfPart, buf);
        if (compare < 0) {
          L = mid + 1;
        } else {
          R = mid;
        }
      }
      segpartptr = partsBySeqId [R];
      if (StringCmp (segpartptr->seqIdOfPart, buf) == 0) {
        if (pos >= segpartptr->from && pos <= segpartptr->to) {
          return segpartptr;
        }
      }
    }
    sip = sip->next;
  }

  return NULL;
}

NLM_EXTERN Int4 LIBCALL SeqMgrMapPartToSegmentedBioseq (BioseqPtr in, Int4 pos, BioseqPtr bsp, SeqIdPtr sip)

{
  BioseqExtraPtr  bspextra;
  SMSeqIdxPtr     currp;
  SMSeqIdxPtr     nextp;
  ObjMgrDataPtr   omdp;
  SMSeqIdxPtr     segpartptr;

  if (in == NULL) return -1;

  /* first check to see if part has been loaded and single map up block installed */

  if (bsp != NULL) {
    omdp = SeqMgrGetOmdpForPointer (bsp);
    if (omdp != NULL) {
      bspextra = (BioseqExtraPtr) omdp->extradata;
      if (bspextra != NULL) {

        /* no need for partsByLoc or partsBySeqId arrays, just use segparthead linked list */

        for (segpartptr = bspextra->segparthead; segpartptr != NULL; segpartptr = segpartptr->next) {
          if (segpartptr->parentBioseq == in) {
            if (pos >= segpartptr->from && pos <= segpartptr->to) {

              /* success, immediate return with mapped up value */

              if (segpartptr->strand == Seq_strand_minus) {
                return segpartptr->cumOffset + (segpartptr->to - pos);
              } else {
                return segpartptr->cumOffset + (pos - segpartptr->from);
              }
            }
          }
        }
      }
    }
  }

  /* otherwise do binary search on segmented bioseq mapping data */

  segpartptr = BinarySearchPartToSegmentMap (in, pos, bsp, sip);
  if (segpartptr == NULL) return -1;

  if (pos >= segpartptr->from && pos <= segpartptr->to) {

    /* install map up block on part, if it has been loaded, to speed up next search */

    if (bsp != NULL) {
      omdp = SeqMgrGetOmdpForPointer (bsp);
      if (omdp != NULL) {
        bspextra = (BioseqExtraPtr) omdp->extradata;
        if (bspextra == NULL) {
          CreateBioseqExtraBlock (omdp, bsp);
          bspextra = (BioseqExtraPtr) omdp->extradata;
        }
        if (bspextra != NULL) {

          /* clean up any old map up info on part */

          for (currp = bspextra->segparthead; currp != NULL; currp = nextp) {
            nextp = currp->next;
            SeqLocFree (currp->slp);
            MemFree (currp->seqIdOfPart);
            MemFree (currp);
          }
          bspextra->segparthead = NULL;
          bspextra->numsegs = 0;
          bspextra->partsByLoc = MemFree (bspextra->partsByLoc);
          bspextra->partsBySeqId = MemFree (bspextra->partsBySeqId);

          /* allocate single map up block */

          currp = MemNew (sizeof (SMSeqIdx));
          if (currp != NULL) {
            currp->slp = AsnIoMemCopy (segpartptr->slp,
                                       (AsnReadFunc) SeqLocAsnRead,
                                       (AsnWriteFunc) SeqLocAsnWrite);
            currp->seqIdOfPart = StringSave (segpartptr->seqIdOfPart);
            currp->parentBioseq = segpartptr->parentBioseq;
            currp->cumOffset = segpartptr->cumOffset;
            currp->from = segpartptr->from;
            currp->to = segpartptr->to;
            currp->strand = segpartptr->strand;
          }

          /* add new map up block to part */

          bspextra->segparthead = currp;
        }
      }
    }

    /* now return offset result */

    if (segpartptr->strand == Seq_strand_minus) {
      return segpartptr->cumOffset + (segpartptr->to - pos);
    } else {
      return segpartptr->cumOffset + (pos - segpartptr->from);
    }
  }
  return -1;
}

