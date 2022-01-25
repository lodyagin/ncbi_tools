/* $Id: ncbilist.hpp,v 1.3 1997/11/12 22:14:33 vakatov Exp $ */
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

 Classes:
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
                            implemented all functions which were N/Impl.

*****************************************************************************/

#ifndef __NCBILIST_HPP__
#define __NCBILIST_HPP__ ncbilist_hpp

#include <ncbibase.hpp>


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* PUBLIC-USE CLASS & FUNCTION DEFINITIONS                                  */
/****************************************************************************/

/**********************************************************************
 * Base classes for double-link list(link, list & iterator)
 **********************************************************************/

class CNcbiBaseLink : public CNcbiObject
{
    friend class CNcbiBaseList;
private: /* these two and Clone() hint is to avoid the chain doubling */
    CNcbiBaseLink(const CNcbiBaseLink&)  { NO_GOOD; }
    CNcbiBaseLink& operator =(const CNcbiBaseLink&)  {
        NO_GOOD;  return *this; }

protected:
    CNcbiBaseLink* m_pPrev;
    CNcbiBaseLink* m_pNext;

public:
    CNcbiBaseLink()  { m_pNext = m_pPrev = NULL; }
    virtual CNcbiObject* Clone() const; /* new Next/Prev = NULL! */
    virtual ~CNcbiBaseLink();

    CNcbiBaseLink* Next() const { return m_pNext; }
    CNcbiBaseLink* Prev() const { return m_pPrev; }
};


class CNcbiBaseList : public CNcbiObject
{
protected:
    CNcbiBaseLink* m_pHead;
    CNcbiBaseLink* m_pTail;
    Uint4 m_iNumItems;

    void Init()  { m_pHead = m_pTail = NULL;  m_iNumItems = 0; }
    CNcbiBaseList& CopyFrom(const CNcbiBaseList& theList);
    ENcbiBoolean AddAfter(CNcbiBaseLink* pAfter, CNcbiBaseLink* pLink);

public:
    /* constructor, assignment, clone & destructor */
    CNcbiBaseList()  { Init(); }
    CNcbiBaseList(CNcbiBaseLink* pLink);
    CNcbiBaseList(const CNcbiBaseList& theList);
    CNcbiBaseList& operator =(const CNcbiBaseList& theList);
    virtual CNcbiObject* Clone() const;
    virtual ~CNcbiBaseList();

    /* member data access functions */
    CNcbiBaseLink* First  () const  { return m_pHead; }
    CNcbiBaseLink* Last   () const  { return m_pTail; }
    Uint4          Size   () const  { return m_iNumItems; }
    ENcbiBoolean   IsEmpty() const  { return ToBoolean(m_iNumItems == 0); }

    /* append new item to the list's head/tail */
    ENcbiBoolean Append (CNcbiBaseLink* pLink);
    ENcbiBoolean Prepend(CNcbiBaseLink* pLink);

    /* exclude item "pLink" from the list */
    CNcbiBaseLink* Remove(CNcbiBaseLink* pLink);

    CNcbiBaseLink* RemoveFirst() { return Remove( First() ); }
    CNcbiBaseLink* RemoveLast () { return Remove( Last () ); }

    /* destroy all items in the list (NOTE:  calls "delete"!)*/
    virtual void Clean();
};


class CNcbiBaseIterator : public CNcbiObject
{
private:
    CNcbiBaseIterator() { m_pList = NULL;  m_pHere = NULL; }

protected:
    CNcbiBaseList* m_pList;
    CNcbiBaseLink* m_pHere;

public:
    /* Constructor, clone & destructor
     * (the default bit-copy copy constructor and assignment are okay)
     */
    CNcbiBaseIterator(const CNcbiBaseList& theList);
    virtual CNcbiObject* Clone() const;
    virtual ~CNcbiBaseIterator();

    /* Member functions */
    CNcbiBaseLink* Next  ();
    CNcbiBaseLink* Cursor() const { return m_pHere; }
    CNcbiBaseLink* Remove();

    ENcbiBoolean   IsNotDone() const { return ToBoolean(m_pHere != NULL); }
    CNcbiBaseList* Container() const { return m_pList; };
    void Reset();
    void Reset(CNcbiBaseList& theList);
};



/**********************************************************************
 * Double-link list, etc. (inherited from the base list classes above)
 **********************************************************************/

class CNcbiListNode : public CNcbiBaseLink
{
private:
    void* m_pData;

    /* these and the Clone() hint is to avoid the chain doubling */
    CNcbiListNode(const CNcbiListNode&) { NO_GOOD; }
    CNcbiListNode& operator =(const CNcbiListNode&) { NO_GOOD; return *this; }
    CNcbiListNode() { NO_GOOD; }

public:
    CNcbiListNode(void* pData) : CNcbiBaseLink() { m_pData = pData; }
    virtual CNcbiObject* Clone() const; /* new Next/Prev = NULL! */
    virtual ~CNcbiListNode();

    CNcbiListNode* Next() const {return (CNcbiListNode*)CNcbiBaseLink::Next();}
    CNcbiListNode* Prev() const {return (CNcbiListNode*)CNcbiBaseLink::Prev();}

    void* Data() const { return m_pData; }
    void  SetData(void *pData) { m_pData = pData; }
};


class CNcbiList : public CNcbiBaseList
{
public:
    CNcbiList() : CNcbiBaseList() {};
    CNcbiList(void* pData);
    CNcbiList(const CNcbiList& theList);
    CNcbiList& operator=(const CNcbiList& theList);

    void* First() const; /* value of the first node of the list */
    void* Last () const; /*           ...last...                */

    ENcbiBoolean Append (void* pData); /* append to the list tail */
    ENcbiBoolean Prepend(void* pData); /*                 ...head */

    /* find & remove first met node with the value equal to "pData"
     * return kNcbiBad if cannot find such a node
     */
    ENcbiBoolean Remove(void* pData);

    void* RemoveFirst();
    void* RemoveLast ();

    static void* Peel(CNcbiListNode* pNode); /* dispose the node; get value */
};


class CNcbiListIterator : public CNcbiBaseIterator
{
public:
    CNcbiListIterator(const CNcbiList& theList);

    void* Next();          /* iterate to the next node;  get its value */
    void* Cursor() const;  /* get curr node value */
    void* Remove();        /* remove curr node from container; get its value */

    CNcbiList* Container() const;
};


class CNcbiStack : public CNcbiList
{
public:
    ENcbiBoolean Push(void*  pData);
    ENcbiBoolean Pop (void*& pData);

    int Depth() const { return Size(); }

    CNcbiList& List() { return (CNcbiList&)*this; }
    const CNcbiList& List() const { return *this; }
};



/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* INLINE FUNCTIONS(IMPLEMENTATION)                                         */
/****************************************************************************/

/* CNcbiBaseList:: member functions
 */

inline CNcbiBaseList::CNcbiBaseList(CNcbiBaseLink* pLink) {
    Init();
    if ( pLink )
        AddAfter(NULL, pLink);
}

inline CNcbiBaseList::CNcbiBaseList(const CNcbiBaseList& theList) {
    Init();
    CopyFrom(theList);
}

inline CNcbiBaseList& CNcbiBaseList::operator =(const CNcbiBaseList& theList) {
    if (&theList == this)
        CopyFrom(theList);
    return *this;
}

inline ENcbiBoolean CNcbiBaseList::Append(CNcbiBaseLink* pLink) {
    return AddAfter(m_pTail, pLink);
}

inline ENcbiBoolean CNcbiBaseList::Prepend(CNcbiBaseLink* pLink) {
    return AddAfter(NULL, pLink);
}


/* CNcbiBaseIterator:: member functions
 */

inline CNcbiBaseIterator::CNcbiBaseIterator(const CNcbiBaseList& theList) {
    m_pList = (CNcbiBaseList*)&theList;
    m_pHere = NULL;
}

inline CNcbiBaseLink* CNcbiBaseIterator::Next() {
    return m_pHere ? (m_pHere = m_pHere->Next()) : 0;
}

inline CNcbiBaseLink* CNcbiBaseIterator::Remove() {
    CNcbiBaseLink* pLink = m_pHere;
    Next();
    return m_pList->Remove(pLink);
}

inline void CNcbiBaseIterator::Reset() {
    if ( m_pList )
        m_pHere = m_pList->First();
}

inline void CNcbiBaseIterator::Reset(CNcbiBaseList& theList) {
    m_pList = &theList;
    Reset();
}


/* CNcbiList:: member functions
 */

inline CNcbiList::CNcbiList(void* pData) {
    if ( pData )
        CNcbiBaseList::Append(new CNcbiListNode(pData));
}

inline CNcbiList::CNcbiList(const CNcbiList& theList)
    : CNcbiBaseList(theList) {};

inline CNcbiList& CNcbiList::operator=(const CNcbiList& theList) {
    return (CNcbiList&)CNcbiBaseList::CopyFrom(theList);
}

inline void* CNcbiList::First() const {
    CNcbiListNode* pNode = (CNcbiListNode*)CNcbiBaseList::First();
    return  pNode ? pNode->Data() : NULL;
}

inline void* CNcbiList::Last() const {
    CNcbiListNode* pNode = (CNcbiListNode*)CNcbiBaseList::Last();
    return  pNode ? pNode->Data() : NULL;
}

inline ENcbiBoolean CNcbiList::Append (void* pData) {
    return CNcbiBaseList::Append(new CNcbiListNode(pData));
}

inline ENcbiBoolean CNcbiList::Prepend(void* pData) {
        return CNcbiBaseList::Prepend(new CNcbiListNode(pData));
}

inline void* CNcbiList::RemoveFirst() {
    return Peel( (CNcbiListNode *)CNcbiBaseList::RemoveFirst() );
}

inline void* CNcbiList::RemoveLast() {
    return Peel( (CNcbiListNode *)CNcbiBaseList::RemoveLast() );
}


/* CNcbiListIterator:: member functions
 */

inline CNcbiListIterator::CNcbiListIterator(const CNcbiList& theList)
    : CNcbiBaseIterator(theList) {}

inline void* CNcbiListIterator::Next() {
    CNcbiListNode* pNode = (CNcbiListNode*)CNcbiBaseIterator::Next();
    return pNode ? pNode->Data() : NULL;
}

inline void* CNcbiListIterator::Cursor() const {
    CNcbiListNode* pNode = (CNcbiListNode*)CNcbiBaseIterator::Cursor();
    return pNode ? pNode->Data() : NULL;
}

inline void* CNcbiListIterator::Remove() {
    return CNcbiList::Peel( (CNcbiListNode*)CNcbiBaseIterator::Remove() );
}

inline CNcbiList* CNcbiListIterator::Container() const {
    return (CNcbiList*)CNcbiBaseIterator::Container();
}


/* CNcbiStack:: member functions
 */

inline ENcbiBoolean CNcbiStack::Push(void*  pData) {
    return Append(pData);
}

inline ENcbiBoolean CNcbiStack::Pop(void*& pData) {
    void *pLast = Last();
    return Remove(pLast) ? (pData = pLast), kNcbiGood : kNcbiBad;
}

#endif /* __NCBILIST_HPP__ */
