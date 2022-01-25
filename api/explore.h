/*  explore.h
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
* File Name:  explore.h
*
* Author:  Jonathan Kans, Jinghui Zhang, James Ostell
*   
* Version Creation Date: 6/30/98
*
* $Revision: 6.15 $
*
* File Description:  Reengineered and optimized exploration functions
*                      to be used for future code
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBI_Explore_
#define _NCBI_Explore_

#ifndef _NCBI_Seqset_
#include <objsset.h>
#endif

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
*
*   SeqMgrBioseqContext, SeqMgrSegmentContext, SeqMgrDescContext and SeqMgrFeatContext
*     are data structures supporting the collection of bioseqs, parts of segmented bioseqs,
*     descriptors, and features, respectively
*
*****************************************************************************/

typedef struct seqmgrbioseqcontext {
  Uint2         entityID;
  Uint2         itemID;
  BioseqPtr     bsp;
  SeqEntryPtr   sep;
  BioseqSetPtr  bssp;
  Pointer       userdata;
                          /* the following fields are for internal use only */
  Pointer       omdp;
  Int2          index;
} SeqMgrBioseqContext, PNTR SeqMgrBioseqContextPtr;

typedef struct seqmgrsegmentcontext {
  Uint2         entityID;
  Uint2         itemID;
  SeqLocPtr     slp;
  BioseqPtr     parent;
  Int4          cumOffset;
  Int4          from;
  Int4          to;
  Uint1         strand;
  Pointer       userdata;
                          /* the following fields are for internal use only */
  Pointer       omdp;
  Int2          index;
} SeqMgrSegmentContext, PNTR SeqMgrSegmentContextPtr;

typedef struct seqmgrdesccontext {
  Uint2         entityID;
  Uint2         itemID;
  ValNodePtr    sdp;
  SeqEntryPtr   sep;
  Uint1         seqdesctype;
  Pointer       userdata;
                          /* the following fields are for internal use only */
  Pointer       omdp;
  Int2          index;
} SeqMgrDescContext, PNTR SeqMgrDescContextPtr;

typedef struct seqmgrfeatcontext {
  Uint2         entityID;
  Uint2         itemID;
  SeqFeatPtr    sfp;
  SeqAnnotPtr   sap;
  CharPtr       label;
  Int2          left;
  Int2          right;
  Uint1         strand;
  Uint2         seqfeattype;
  Uint2         featdeftype;
  Int2          numivals;
  Int4Ptr       ivals;
  Pointer       userdata;
                          /* the following fields are for internal use only */
  Pointer       omdp;
  Int2          index;
} SeqMgrFeatContext, PNTR SeqMgrFeatContextPtr;

/*****************************************************************************
*
*   SeqMgrIndexFeatures builds indices of sorted features for all bioseqs in an
*     entity, makes explicit connections from a protein bioseq to its best protein
*     feature and to the CDS feature encoding it, can be called given an entityID
*     or a BioseqPtr or SeqEntryPtr, and returns the entityID
*
*****************************************************************************/

NLM_EXTERN Uint2 LIBCALL SeqMgrIndexFeatures PROTO((Uint2 entityID, Pointer ptr, Uint1 labeltype));

/*****************************************************************************
*
*   SeqMgrGetBestProteinFeature and SeqMgrGetCDSgivenProduct take a protein
*     bioseq to get the best protein feature or encoding CDS
*   SeqMgrGetRNAgivenProduct takes an mRNA (cDNA) bioseq and gets the encoding
*     mRNA feature on the genomic bioseq
*
*****************************************************************************/

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetBestProteinFeature PROTO((BioseqPtr bsp,
                                                                 SeqMgrFeatContext PNTR context));

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetCDSgivenProduct PROTO((BioseqPtr bsp,
                                                              SeqMgrFeatContext PNTR context));

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetRNAgivenProduct PROTO((BioseqPtr bsp,
                                                              SeqMgrFeatContext PNTR context));

/*****************************************************************************
*
*   To find the best gene feature, first call SeqMgrGetGeneXref, and if it is not
*     NULL call SeqMgrGeneIsSuppressed, otherwise call SeqMgrGetOverlappingGene
*   If desired, place a SeqMgrFeatContext data structure on the stack, and pass
*     in &context as the second parameter to SeqMgrGetOverlappingGene
*
*****************************************************************************/

NLM_EXTERN GeneRefPtr LIBCALL SeqMgrGetGeneXref PROTO((SeqFeatPtr sfp));

NLM_EXTERN Boolean LIBCALL SeqMgrGeneIsSuppressed PROTO((GeneRefPtr grp));

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetOverlappingGene PROTO((SeqFeatPtr sfp,
                                                              SeqMgrFeatContext PNTR context));

/*****************************************************************************
*
*   SeqMgrGetOverlappingPub returns the overlapping publication feature
*   If desired, place a SeqMgrFeatContext data structure on the stack, and pass
*     in &context as the second parameter to SeqMgrGetOverlappingPub
*
*****************************************************************************/

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetOverlappingPub PROTO((SeqFeatPtr sfp,
                                                             SeqMgrFeatContext PNTR context));

/*****************************************************************************
*
*   SeqMgrGetOverlappingSource returns the overlapping biosource feature
*   If desired, place a SeqMgrFeatContext data structure on the stack, and pass
*     in &context as the second parameter to SeqMgrGetOverlappingSource
*
*****************************************************************************/

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetOverlappingSource PROTO((SeqFeatPtr sfp,
                                                                SeqMgrFeatContext PNTR context));

/*****************************************************************************
*
*   Replacements for BioseqContext functions using bioseq feature indices, returning
*     the next (sorted) feature or (bioseq to highest set) descriptor pointer
*   The SeqMgrDescContext or SeqMgrFeatContext data structures should be on the
*     calling function's stack, and are passed as &context to the context function.
*   Passing NULL for curr in the first call initializes the context structure, and
*     the functions return NULL at the end of the list
*   If the choice parameters are 0, every feature or descriptor is returned
*   It is expected that these calls would be flanked by BioseqLock and BioseqUnlock,
*     so object manager reload could ensure that pointers are valid within the loop,
*     since the pointers are what drive these functions
*   The Explore functions below offer more flexibility than these Context functions
*
*****************************************************************************/

NLM_EXTERN ValNodePtr LIBCALL SeqMgrGetNextDescriptor PROTO((BioseqPtr bsp, ValNodePtr curr,
                                                             Uint1 seqDescChoice,
                                                             SeqMgrDescContext PNTR context));

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetNextFeature PROTO((BioseqPtr bsp, SeqFeatPtr curr,
                                                          Uint1 seqFeatChoice, Uint1 featDefChoice,
                                                          SeqMgrFeatContext PNTR context));

/*****************************************************************************
*
*   Callback types for SeqMgrExploreBioseqs, SeqMgrExploreSegments,
*     SeqMgrExploreDescriptors, and SeqMgrExploreFeatures
*
*****************************************************************************/

typedef Boolean (LIBCALLBACK *SeqMgrBioseqExploreProc) PROTO((BioseqPtr bsp, SeqMgrBioseqContextPtr context));

typedef Boolean (LIBCALLBACK *SeqMgrSegmentExploreProc) PROTO((SeqLocPtr slp, SeqMgrSegmentContextPtr context));

typedef Boolean (LIBCALLBACK *SeqMgrDescExploreProc) PROTO((ValNodePtr sdp, SeqMgrDescContextPtr context));

typedef Boolean (LIBCALLBACK *SeqMgrFeatExploreProc) PROTO((SeqFeatPtr sfp, SeqMgrFeatContextPtr context));

/*****************************************************************************
*
*   SeqMgrExploreBioseqs, SeqMgrExploreSegments, SeqMgrExploreDescriptors, and
*     SeqMgrExploreFeatures use the bioseq feature indices to quickly present
*     desired items to the user-supplied callback function, stopping if the callback
*     returns FALSE
*   In contrast to the SeqMgrGetNext functions, the SeqMgrExplore function callbacks pass
*     a pointer to the SeqMgr[Bioseq/Segment/Desc/Feat]Context data structures held by the
*     explore function, not on the calling function's stack
*   If the filter parameters are NULL, every feature or descriptor is returned, otherwise
*     the array lengths should be SEQDESCR_MAX, SEQFEAT_MAX, and FEATDEF_MAX, and the
*     elements are from the Seq_descr_, SEQFEAT_, and FEATDEF_ lists
*   It is expected that these calls would be flanked by BioseqLock and BioseqUnlock,
*     so object manager reload could ensure that pointers are valid within the loop,
*     but these explore functions can work on cached-out records
*
*****************************************************************************/

NLM_EXTERN Boolean LIBCALL SeqMgrExploreBioseqs PROTO((Uint2 entityID, Pointer ptr, Pointer userdata,
                                                       SeqMgrBioseqExploreProc userfunc,
                                                       Boolean nucs, Boolean prots, Boolean parts));

NLM_EXTERN Boolean LIBCALL SeqMgrExploreSegments PROTO((BioseqPtr bsp, Pointer userdata,
                                                        SeqMgrSegmentExploreProc userfunc));

NLM_EXTERN Boolean LIBCALL SeqMgrExploreDescriptors PROTO((BioseqPtr bsp, Pointer userdata,
                                                          SeqMgrDescExploreProc userfunc,
                                                          BoolPtr seqDescFilter));

NLM_EXTERN Boolean LIBCALL SeqMgrExploreFeatures PROTO((BioseqPtr bsp, Pointer userdata,
                                                        SeqMgrFeatExploreProc userfunc,
                                                        SeqLocPtr locationFilter,
                                                        BoolPtr seqFeatFilter, BoolPtr featDefFilter));

/*****************************************************************************
*
*   SeqMgrGetDesiredDescriptor and SeqMgrGetDesiredFeature return a descriptor
*     or feature given either an itemID, a position index, or the feature or
*     descriptor pointer itself, using whichever parameter is not 0 (or NULL)
*   In order to obtain index information associated with the desired descriptor
*     or feature, place a SeqMgrDescContext or SeqMgrFeatContext data structure
*     on the stack, and pass in &context as the last parameter
*
*****************************************************************************/

NLM_EXTERN ValNodePtr LIBCALL SeqMgrGetDesiredDescriptor PROTO((BioseqPtr bsp, Uint2 itemID,
                                                                Int2 index, ValNodePtr sdp,
                                                                SeqMgrDescContext PNTR context));

NLM_EXTERN SeqFeatPtr LIBCALL SeqMgrGetDesiredFeature PROTO((BioseqPtr bsp, Uint2 itemID,
                                                             Int2 index, SeqFeatPtr sfp,
                                                             SeqMgrFeatContext PNTR context));

/* the following functions are not frequently called by applications */

/*****************************************************************************
*
*   SeqMgrFeaturesAreIndexed returns the last time feature indices were built,
*     with 0 meaning that indices are not present on the entity
*
*****************************************************************************/

NLM_EXTERN time_t LIBCALL SeqMgrFeaturesAreIndexed PROTO((Uint2 entityID));

/*****************************************************************************
*
*   SeqMgrClearFeatureIndexes clears feature indices for an entity given an
*     entityID or a BioseqPtr or SeqEntryPtr
*
*****************************************************************************/

NLM_EXTERN Boolean LIBCALL SeqMgrClearFeatureIndexes PROTO((Uint2 entityID, Pointer ptr));



#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _NCBI_Explore_ */

