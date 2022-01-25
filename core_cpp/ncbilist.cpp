/* $Id: ncbilist.cpp,v 1.2 1997/10/27 18:57:18 vakatov Exp $ */
/*****************************************************************************

 Description:  NCBI basic lists & iterators

 Authors:  Grigoriy Starchenko, Denis Vakatov

***************************************************************************

                       PUBLIC DOMAIN NOTICE
           National Center for Biotechnology Information

 This software/database is a "United States Government Work" under the
 terms of the United States Copyright Act.  It was written as part of    
 the author's official duties as a United States Government employee
 and thus cannot be copyrighted.  This software/database is freely
 available to the public for use. The National Library of Medicine and
 the U.S. Government have not placed any restriction on its use or
 reproduction.

 Although all reasonable efforts have been taken to ensure the accuracy
 and reliability of the software and data, the NLM and the U.S.
 Government do not and cannot warrant the performance or results that
 may be obtained by using this software or data. The NLM and the U.S.
 Government disclaim all warranties, express or implied, including
 warranties of performance, merchantability or fitness for any
 particular purpose.

 Please cite the author in any work or product based on this material.

***************************************************************************

 Implementation of Classes:
   "CNcbiBaseLink     : CNcbiObject"
   "CNcbiBaseList     : CNcbiObject"
   "CNcbiBaseIterator : CNcbiObject"

   "CNcbiListNode     : CNcbiBaseLink"
   "CNcbiList         : CNcbiBaseList"
   "CNcbiListIterator : CNcbiBaseIterator"

   "CNcbiStack        : CNcbiList"


 Modification History:
   Sep  9 1997 (grisha)  -- originally written
   Oct 27 1997 (vakatov) -- moved more code to inline functions
                            implemented all functions which were not impl.

*****************************************************************************/

#include <ncbilist.hpp>


/********************************************************************/
/* CNcbiBaseLink implementation */
/********************************************************************/

CNcbiObject*
/*FCN*/CNcbiBaseLink::Clone() const
{
    return (CNcbiObject*)(new CNcbiBaseLink());
}

/*FCN*/CNcbiBaseLink::~CNcbiBaseLink()
{
    m_pNext = NULL;
    m_pPrev = NULL;
}


/********************************************************************/
/* CNcbiBaseList implementation */
/********************************************************************/
CNcbiBaseList&
/*FCN*/CNcbiBaseList::CopyFrom(const CNcbiBaseList& theList)
{
    Clean();
    for (CNcbiBaseLink *pLink = theList.First();
         pLink;  pLink = pLink->Next()) {
        Append( (CNcbiBaseLink *)pLink->Clone() );
    }
    return *this;
}

CNcbiObject*
/*FCN*/CNcbiBaseList::Clone() const
{
    return (CNcbiObject*)(new CNcbiBaseList(*this));
}


/*FCN*/CNcbiBaseList::~CNcbiBaseList()
{
    Clean();
}


CNcbiBaseLink*
/*FCN*/CNcbiBaseList::Remove(CNcbiBaseLink* pLink)
{
    if ( !pLink )
        return NULL;

    CNcbiBaseLink *pBefore = pLink->m_pPrev;
    CNcbiBaseLink *pAfter  = pLink->m_pNext;
    if ( pBefore ) {  /* any elements before? */
        pBefore->m_pNext = pLink->m_pNext;
    } else {               /* this is a head of the list */
        m_pHead = pAfter;  /* make Head correction */
    }
    if ( pAfter ) {   /* any elements after? */
        pAfter->m_pPrev = pBefore;
    } else {
        m_pTail = pBefore;   /* make Tail correction */
    }
    pLink->m_pPrev = pLink->m_pNext = NULL; 

    m_iNumItems--;
    return pLink;
}


ENcbiBoolean
/*FCN*/CNcbiBaseList::AddAfter (
  CNcbiBaseLink* pAfter,
  CNcbiBaseLink* pLink
){
    if ( !pLink )
        return kNcbiBad;

    if ( pAfter ) {    /* insert new item */
        CNcbiBaseLink *pBefore;
        pLink->m_pPrev = pAfter;
        pLink->m_pNext = pBefore = pAfter->m_pNext;
        pAfter->m_pNext = pLink;
        if ( pBefore ) {
            pBefore->m_pPrev = pLink;
        }
    } else {           /* add item to the head */
        pLink->m_pPrev = NULL;
        pLink->m_pNext = m_pHead;
        if ( m_pHead ) {
            m_pHead->m_pPrev = pLink;
        }
        m_pHead = pLink;
    }
    if (m_pTail == pAfter) { /* check tail for correction */
        m_pTail = pLink;
    }

    m_iNumItems++;
    return kNcbiGood;
}

void
/*FCN*/CNcbiBaseList::Clean()
{
    for (CNcbiBaseLink *pLink = First();  pLink;  pLink = First()) {
        delete Remove( pLink );
    }
    ASSERT ( !m_pHead && !m_pTail && !m_iNumItems );
}


/********************************************************************/
/* CNcbiBaseIterator implementation */
/********************************************************************/

/*FCN*/CNcbiBaseIterator::~CNcbiBaseIterator()
{
    m_pList = NULL;
    m_pHere = NULL;
}

CNcbiObject*
/*FCN*/CNcbiBaseIterator::Clone() const {
    CNcbiBaseIterator *iter = new CNcbiBaseIterator;
    if ( iter )
        *iter = *this;
    return iter;
}


/********************************************************************/
/* CNcbiListNode implementation */
/********************************************************************/

CNcbiObject*
/*FCN*/CNcbiListNode::Clone() const
{
    return new CNcbiListNode( this->Data() );
}

/*FCN*/CNcbiListNode::~CNcbiListNode()
{
    m_pData = NULL;
}


/********************************************************************/
/* CNcbiList implementation */
/********************************************************************/

void*
/*FCN*/CNcbiList::Peel(CNcbiListNode* pNode)
{
    if ( !pNode )
        return NULL;
    ASSERT ( !pNode->Prev()  &&  !pNode->Next() );

    void* pData = pNode->Data();
    delete pNode;
    return pData;
}

ENcbiBoolean
/*FCN*/CNcbiList::Remove(void* pData)
{
    for (CNcbiBaseLink* pLink = CNcbiBaseList::First();
         pLink;  pLink = pLink->Next()) {
        if (((CNcbiListNode*)pLink)->Data() == pData) {
            Peel( (CNcbiListNode*)CNcbiBaseList::Remove(pLink) );
            return kNcbiGood;
        }
    }

    return kNcbiBad;
}

/* EOF */
