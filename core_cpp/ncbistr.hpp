/* $Id: ncbistr.hpp,v 1.11 1997/12/02 21:50:28 vakatov Exp $
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
* Authors:  Grigoriy Starchenko, Denis Vakatov
*
* File Description:
*    "CNcbiString : CNcbiObject"
*       -- NCBI string class and base NCBI string class
*       -- (Note:  to provide a subset of ANSI string class interface)
*
*    "CNcbiStringList : CNcbiList"
*    "CNcbiStringListIterator : CNcbiListIterator"
*       -- container & iterator classes for the NCBI string class
*
*    "CNcbiBaseStr":
*       -- base NCBI string class w/reference counter
*
* Modifications:
* --------------------------------------------------------------------------
* Sep  9 1997 (grisha)  -- originally written
* Oct 27 1997 (vakatov) -- redesigned:
*                             added CNcbiBaseStr
*                             reimplemented CNcbiString
*                             moved more code to inline functions
*
* ==========================================================================
*/

#ifndef __NCBISTR_HPP__
#define __NCBISTR_HPP__

#include <ncbilist.hpp>

// by default, do not reserve extra space for the string character buffer
// (default: #define DEF_SIZE_POLICY eReserveSize)
// #define DEF_SIZE_POLICY eExactSize


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* PUBLIC-USE CLASS & FUNCTION DEFINITIONS                                  */
/****************************************************************************/

/****************************************************************************
 * CLASS CNcbiString
 ****************************************************************************/

class CNcbiBaseStr;

class CNcbiString : public CNcbiObject
{
private:
    CNcbiBaseStr *m_BS;

    void GetReadyForUpdate(Uint4 length=0);
public:
    /* Constructors */
    CNcbiString(const CNcbiString& ns);
    CNcbiString(const Char *cs=NULL, Uint4 length=0);

    /* Clone(see also CNcbiObject) */
    virtual CNcbiObject *Clone() const;

    /* Destructor(see also CNcbiObject) */
    virtual ~CNcbiString();

    /* Properties */
    Uint4        BufferSize () const;  // size of the allocated string buffer
    Uint4        Length     () const;  // length of the string
    ENcbiBoolean Empty      () const;  // empty string ("") or NULL string
    ENcbiBoolean Null       () const;  // NULL string

    /* Operator assign */
    CNcbiString& operator=(const CNcbiString& ns);
    CNcbiString& operator=(const Char* cs);  // (makes its own copy of "cs")
    CNcbiString& operator=(Char ch); 

    /* Access */
    Char operator[](Uint4 pos) const;

    /* Simple edit */
    CNcbiString& SetChar   (Uint4 pos, Char ch);
    CNcbiString& DeleteChar(Uint4 pos, Uint4 n_del=1);
    CNcbiString& InsertChar(Uint4 pos, Char ch, Uint4 n_ins=1);
    CNcbiString& InsertStr (Uint4 pos, const Char *cs, Uint4 n_ins=UINT4_MAX);

    /* Direct data access -- be very careful here!
     */
    // the pointer returned by the following 3 funcs is guaranteed to
    // be valid only until *any* change is made in the string; after any
    // change, the pointer can become invalid!
    operator const Char*() const; // these two return pointer to the
    const Char*    Data () const; //  wrapped C string(can return NULL)
    const Char*    C_Str() const; // (return empty C string in the "NULL case")
  
    // get the string character buffer for update:
    // dont write out of "BufferSize()";  don't use any(but the "BufferSize()")
    // of the CNcbiString functions until you stop your direct
    // changes and call "CommitBuffer()"
    Char *FreezeBuffer();
    // put (updated) buffer back and make it in-sync with other
    // class members;  if "length" != UINT4_MAX specified then
    // use only first "length" characters
    void  CommitBuffer(Uint4 length=UINT4_MAX);

    /* Member functions */
    CNcbiString& Append(const CNcbiString& ns);
    CNcbiString& Append(const Char *cs, Uint4 len=UINT4_MAX);
    CNcbiString& Append(Char ch, Uint4 len=1);

    CNcbiString& operator+=(const CNcbiString& ns);
    CNcbiString& operator+=(const Char *cs);
    CNcbiString& operator+=(Char ch);

    void Swap(CNcbiString& ns);

    /* safe variant of "strlen()" -- do not scan for '\0' after "max_len"
     * position;  return pos. <= "max_len" */
    static Uint4 strnlen(const Char *buf, Uint4 max_len);

    /* return 0 if first "n" symbols(up to the '\0' symbol) of the two
     * strings("s1" and "s2") are equal(or if both of them are NULL);
     * return 1-based position of the first non-matching character;
     * the returned position is negative if s2 is "bigger" than s1;
     * NULL string is less than any non-NULL string(ret. {-1, +1});
     */
    static Int4 Compare(const Char *s1, const Char *s2, Int4 n=INT4_MAX);
    Int4 Compare(const Char *cs, Int4 n=INT4_MAX) const;

    CNcbiString& Clean();  // set content to NULL
    CNcbiString& ToUpper();
    CNcbiString& ToLower();
    CNcbiString UrlVersion () const;
    CNcbiString HtmlVersion() const;

#ifdef NCBISTR_TODO
    ENcbiBoolean ToInt(Int4& value) const;
    ENcbiBoolean StartsWith(const Char *cs) const;
    CNcbiString  operator()(Uint4 pos, Uint4 len=UINT4_MAX) const;
    CNcbiString  Field(Uint4 start_pos, const Char* csep) const;

    /* search forward/backward for the substring "cs" or character "ch";
     * return the 0-based position of found substring or character;
     * return NCBISTR_ERR if cannot find any */
#define NCBISTR_ERR ((Uint4)~0)
    Uint4 Index(const Char *cs, Uint4 start_pos=0) const;
    Uint4 Index(Char ch, Uint4 start_pos=0) const;
    Uint4 IndexBack(const Char *cs, Uint4 start_pos=UINT4_MAX) const;
    Uint4 IndexBack(Char ch, Uint4 start_pos=UINT4_MAX) const;

    /* starting from position "start_pos", replace no more than "max_replace"
     * occurences of substring "search" by string "replace" */
    Uint4 Replace(const Char *search, const Char *replace,
                  Uint4 start_pos=0, Uint4 max_replace=UINT4_MAX);
#endif /* NCBISTR_TODO */

    /* trim leading and/or trailing symbols specified by "charset" */
    static const Char *kDefaultTrimCharset; // " \t\n\r"
    CNcbiString& Trim     (const Char *charset=kDefaultTrimCharset);
    CNcbiString& TrimRight(const Char *charset=kDefaultTrimCharset);
    CNcbiString& TrimLeft (const Char *charset=kDefaultTrimCharset);
    /* remove all symbols enlisted in "charset" */
    CNcbiString& Prune(const Char *charset=kDefaultTrimCharset);
};

/* binary operators over the "CNcbiString" and "const Char*" */
inline CNcbiString  operator +(const CNcbiString& ns,  const Char *cs);
inline CNcbiString  operator +(const Char *cs,         const CNcbiString& ns);
inline CNcbiString  operator +(Char  ch,               const CNcbiString& ns);
inline CNcbiString  operator +(const CNcbiString& ns,  Char ch);
inline CNcbiString  operator +(const CNcbiString& ns1, const CNcbiString& ns2);

inline ENcbiBoolean operator==(const CNcbiString& ns,  const Char *cs);
inline ENcbiBoolean operator==(const Char *cs,         const CNcbiString& ns);
inline ENcbiBoolean operator==(const CNcbiString& ns1, const CNcbiString& ns2);
inline ENcbiBoolean operator!=(const CNcbiString& ns,  const Char *cs);
inline ENcbiBoolean operator!=(const Char *cs,         const CNcbiString& ns);
inline ENcbiBoolean operator!=(const CNcbiString& ns1, const CNcbiString& ns2);
inline ENcbiBoolean operator <(const CNcbiString& ns,  const Char *cs);
inline ENcbiBoolean operator <(const Char *cs,         const CNcbiString& ns);
inline ENcbiBoolean operator <(const CNcbiString& ns1, const CNcbiString& ns2);
inline ENcbiBoolean operator<=(const CNcbiString& ns,  const Char *cs);
inline ENcbiBoolean operator<=(const Char *cs,         const CNcbiString& ns);
inline ENcbiBoolean operator<=(const CNcbiString& ns1, const CNcbiString& ns2);
inline ENcbiBoolean operator >(const CNcbiString& ns,  const Char *cs);
inline ENcbiBoolean operator >(const Char *cs,         const CNcbiString& ns);
inline ENcbiBoolean operator >(const CNcbiString& ns1, const CNcbiString& ns2);
inline ENcbiBoolean operator>=(const CNcbiString& ns,  const Char *cs);
inline ENcbiBoolean operator>=(const Char *cs,         const CNcbiString& ns);
inline ENcbiBoolean operator>=(const CNcbiString& ns1, const CNcbiString& ns2);



/****************************************************************************
 * CLASSes CNcbiStringList & CNcbiStringListIterator
 * Container & iterator for CNcbiString objects
 ****************************************************************************/

class CNcbiStringList : public CNcbiList
{
public:
    CNcbiStringList() : CNcbiList() {};
    ~CNcbiStringList();

    /* head/tail */
    CNcbiString *First() const;
    CNcbiString *Last () const;

    /* Add string to the head or tail of the list
     * Return NULL on failure
     */
    CNcbiString *Append(const CNcbiString& ns);
    CNcbiString *Append(const Char *cs);
    CNcbiString *Prepend(const CNcbiString& ns);
    CNcbiString *Prepend(const Char *cs);

    /* Unlink "ns" from the list
     * (NOTE: "ns" must be a value returned by Append/Prepend functions)
     * The deallocation of returned pointer is up to caller
     * Return  kNcbiBad if cannot find one
     */
    ENcbiBoolean Remove(const CNcbiString *ns);
    CNcbiString *RemoveFirst();
    CNcbiString *RemoveLast ();

    /* Unlink "ns" from the list, then destroy(deallocate) "ns"
     * Return  kNcbiBad if cannot find one
     */
    ENcbiBoolean Destroy(CNcbiString*& ns);
    ENcbiBoolean DestroyFirst();
    ENcbiBoolean DestroyLast ();

    /* destroy all items in the list */
    virtual void Clean();
};


class CNcbiStringListIterator : public CNcbiListIterator
{
public:
    CNcbiStringListIterator(const CNcbiStringList& theList)
        : CNcbiListIterator(theList) {};

    CNcbiString* Next() {
        return (CNcbiString*)CNcbiListIterator::Next  (); }
    CNcbiString* Cursor() const {
        return (CNcbiString*)CNcbiListIterator::Cursor(); }
    operator CNcbiString *() const { return Cursor(); }
    CNcbiString* Remove() {
        return (CNcbiString*)CNcbiListIterator::Remove(); }
    ENcbiBoolean Destroy() {
        CNcbiString* ns = Remove();
        return ns ? delete ns, kNcbiTrue : kNcbiFalse;
    } 

    CNcbiStringList* Container() const {
        return (CNcbiStringList*)CNcbiListIterator::Container(); };
};



/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* INTERNAL-USE CLASS & FUNCTION DEFINITIONS                                */
/****************************************************************************/

/****************************************************************************
 * CLASS CNcbiBaseStr
 * Wraps a standard '\0'-terminated C single-byte string(char *);  adds:
 *  reference counter,
 *  string length,
 *  size of allocated buffer;
 *  set of basic methods -- create, clone, dispose, etc.
 ****************************************************************************/

class CNcbiBaseStr
{
public:
#ifndef DEF_SIZE_POLICY
#define DEF_SIZE_POLICY eReserveSize
#endif
    enum ESizePolicy {
        eExactSize,
        eReserveSize
    };

private:
    struct SSmallAttr { // to be used if the string is small(<=255 chars)
        Uint2 ref_count;
        Uint1 size;
        Uint1 length;
    };

    struct SBigAttr { // to be used if the string is large(>255 chars)
        Uint4 ref_count;
        Uint4 size;
        Uint4 length;
        Uint4 is_small;  // use the above fields if "is_small" = 0
    };

    Char m_Str[1]; // the character string itself

    Uint4 IsSmall() const { return *((Uint4 *)this - 1); }

    SBigAttr& BigAttr() {
        ASSERT ( !IsSmall() );  return *((SBigAttr *)this - 1); }
    const SBigAttr& BigAttr() const {
        ASSERT ( !IsSmall() );  return *((SBigAttr *)this - 1); }

    SSmallAttr& SmallAttr() {
        ASSERT ( IsSmall() );  return *((SSmallAttr *)this - 1); }
    const SSmallAttr& SmallAttr() const {
        ASSERT ( IsSmall() );  return *((SSmallAttr *)this - 1); }

    void *operator new(size_t,
                       const Char *str=0,
                       Uint4 str_length=0,
                       ESizePolicy sz_policy=DEF_SIZE_POLICY);
    void operator delete(void *ptr);

    // special case -- corresponds to "(char *)NULL" in assignment and
    // comparison;  otherwise behaves just like an empty string
    friend struct SFakeNullStr;
    static CNcbiBaseStr *sm_NullStr;

    static Uint4 sm_Counter;      // counts all CNcbiBaseStr class instances
    static Uint4 sm_SmallCounter; // ... small, alloc by "new" ...
    static Uint4 sm_BigCounter;   // ... big,   alloc by "new" ...

    // prohibited for implicit use by outside users
    CNcbiBaseStr() { sm_Counter++; }
    CNcbiBaseStr(const CNcbiBaseStr&) { NO_GOOD;  sm_Counter++; };
    CNcbiBaseStr& operator =(const CNcbiBaseStr&) {
        NO_GOOD;  return *this; }
    ~CNcbiBaseStr() {};

public:
    static CNcbiBaseStr *Create(const Char *str=NULL,
                                Uint4 length=0,
                                ESizePolicy sz_policy=DEF_SIZE_POLICY) {
        return new(str, length, sz_policy) CNcbiBaseStr;
    }
    CNcbiBaseStr *Clone();
    CNcbiBaseStr *Realloc(Uint4 length, ESizePolicy sz_policy=DEF_SIZE_POLICY);
    void Release() { delete this; }

    Uint4 RefCount() const {
        return IsSmall() ? SmallAttr().ref_count : BigAttr().ref_count; }
    Uint4 BufferSize() const {
        return IsSmall() ? SmallAttr().size : BigAttr().size; }
    Uint4 Length() const {
        return IsSmall() ? SmallAttr().length : BigAttr().length; }
    void SetLength(Uint4 length) {
        ASSERT ( length < BufferSize() );
        m_Str[length] = '\0';
        ASSERT ( length == strlen(m_Str) );
        if ( IsSmall() )  SmallAttr().length = (Uint1)length;
        else BigAttr().length = length; }

    Char  operator [](Uint4 pos) const {
        ASSERT ( pos <= Length() );  return m_Str[pos]; }
    Char& operator [](Uint4 pos) {
        ASSERT ( pos <= Length() );  return m_Str[pos]; }

    const Char* Str() const { return m_Str; }
    Char* Str() { return m_Str; }

    ENcbiBoolean Empty() const {
        return ToBoolean( !m_Str[0] ); }
    ENcbiBoolean Null() const {
        return ToBoolean(this == sm_NullStr); }

    static Uint4 nInstances() { return sm_Counter; } //for debugging purporses 
};


/****************************************************************************
 * INLINE FUNCTIONS(IMPLEMENTATION)
 ****************************************************************************/

/* CNcbiString:: member functions
 */

inline CNcbiString::CNcbiString(const CNcbiString& ns) {
    m_BS = ns.m_BS->Clone();
}

inline CNcbiString::CNcbiString(const Char *cs, Uint4 length) {
    m_BS = CNcbiBaseStr::Create(cs, length);
}

inline Uint4 CNcbiString::BufferSize() const {
    return m_BS->BufferSize();
}

inline Uint4        CNcbiString::Length() const { return m_BS->Length(); }
inline ENcbiBoolean CNcbiString::Empty () const { return m_BS->Empty (); }
inline ENcbiBoolean CNcbiString::Null  () const { return m_BS->Null  (); }
  
inline CNcbiString& CNcbiString::operator=(const CNcbiString& ns) {
    if (m_BS != ns.m_BS) {
        m_BS->Release();
        m_BS = ns.m_BS->Clone();
    }
    return *this;
}

inline CNcbiString& CNcbiString::operator=(const Char* cs) {
    m_BS->Release();
    m_BS = CNcbiBaseStr::Create( cs );
    return *this;
}

inline CNcbiString& CNcbiString::operator=(Char ch) {
    m_BS->Release();
    m_BS = CNcbiBaseStr::Create(&ch, 1);
    return *this;
}

inline Char CNcbiString::operator[](Uint4 pos) const { return (*m_BS)[pos]; }

inline CNcbiString::operator const Char*() const {
    return Null() ? 0 : m_BS->Str();
}
inline const Char* CNcbiString::Data()  const {
    return Null() ? 0 : m_BS->Str();
}
inline const Char* CNcbiString::C_Str() const {
    return m_BS->Str();
}

inline CNcbiString&  CNcbiString::SetChar(Uint4 pos, Char ch) {
    GetReadyForUpdate();
    (*m_BS)[pos] = ch;
    return *this;
}

inline Char *CNcbiString::FreezeBuffer() {
    GetReadyForUpdate();
    return m_BS->Str();
}

inline void CNcbiString::CommitBuffer(Uint4 length) {
    ASSERT ( m_BS->RefCount() == 1 );
    ASSERT ( length == UINT4_MAX  ||  length < BufferSize() );
    m_BS->SetLength( strnlen(C_Str(), length) );
}

inline CNcbiString& CNcbiString::operator+=(const CNcbiString& ns) {
    return Append( ns );
}
inline CNcbiString& CNcbiString::operator+=(const Char *cs) {
    return Append( cs );
}
inline CNcbiString& CNcbiString::operator+=(Char ch) {
    return Append( ch );
}

inline void CNcbiString::Swap(CNcbiString& ns) {
    CNcbiBaseStr *nbstr = m_BS;
    m_BS = ns.m_BS;
    ns.m_BS = nbstr;
}

inline  CNcbiString& CNcbiString::Clean() {
    return (*this = (const Char*)0);
}

inline Int4 CNcbiString::Compare(const Char *cs, Int4 n) const {
    return CNcbiString::Compare((const Char*)(*this), cs, n);
}


/* CNcbiString external inline functions(binary operators)
 */

inline CNcbiString operator+(const CNcbiString& ns, const Char *cs) {
    if (!cs || !*cs)
        return ns;
    return CNcbiString(ns, ns.Length() + strlen(cs)).Append(cs);
}
inline CNcbiString operator+(const Char *cs, const CNcbiString& ns) {
    if (!cs || !*cs)
        return ns;
    return CNcbiString(cs, strlen(cs) + ns.Length()).Append(ns);
}
inline CNcbiString operator+(Char ch, const CNcbiString& ns) {
    Char cs[2];  cs[0] = ch;  cs[1] = '\0';
    return cs + ns;
}
inline CNcbiString operator+(const CNcbiString& ns, Char ch) {
    return CNcbiString(ns, ns.Length()+1).Append(ch);
}
inline CNcbiString operator+(const CNcbiString& ns1, const CNcbiString& ns2) {
    return CNcbiString(ns1, ns1.Length() + ns2.Length()).Append(ns2);
}

inline ENcbiBoolean operator==(const CNcbiString& ns, const Char *cs) {
    return ToBoolean(CNcbiString::Compare(ns, cs) == 0);
}
inline ENcbiBoolean operator==(const Char *cs, const CNcbiString& ns) {
    return ToBoolean(CNcbiString::Compare(cs, ns) == 0);
}
inline ENcbiBoolean operator==(const CNcbiString& ns1, const CNcbiString& ns2){
    return ToBoolean(CNcbiString::Compare(ns1, (const Char*)ns2) == 0);
}

inline ENcbiBoolean operator!=(const CNcbiString& ns, const Char *cs) {
    return ToBoolean(CNcbiString::Compare(ns, cs) != 0);
}
inline ENcbiBoolean operator!=(const Char *cs, const CNcbiString& ns) {
    return ToBoolean(CNcbiString::Compare(cs, ns) != 0);
}
inline ENcbiBoolean operator!=(const CNcbiString& ns1, const CNcbiString& ns2){
    return ToBoolean(CNcbiString::Compare(ns1, (const Char*)ns2) != 0);
}

inline ENcbiBoolean operator<(const CNcbiString& ns, const Char *cs) {
    return ToBoolean(CNcbiString::Compare(ns, cs) < 0);
}
inline ENcbiBoolean operator<(const Char *cs, const CNcbiString& ns) {
    return ToBoolean(CNcbiString::Compare(cs, ns) < 0);
}
inline ENcbiBoolean operator<(const CNcbiString& ns1, const CNcbiString& ns2) {
    return ToBoolean(CNcbiString::Compare(ns1, (const Char*)ns2) < 0);
}

inline ENcbiBoolean operator<=(const CNcbiString& ns, const Char *cs) {
    return ToBoolean(CNcbiString::Compare(ns, cs) <= 0);
}
inline ENcbiBoolean operator<=(const Char *cs, const CNcbiString& ns) {
    return ToBoolean(CNcbiString::Compare(cs, ns) <= 0);
}
inline ENcbiBoolean operator<=(const CNcbiString& ns1, const CNcbiString& ns2){
    return ToBoolean(CNcbiString::Compare(ns1, (const Char*)ns2) <= 0);
}

inline ENcbiBoolean operator>(const CNcbiString& ns, const Char *cs) {
    return ToBoolean(CNcbiString::Compare(ns, cs) > 0);
}
inline ENcbiBoolean operator>(const Char *cs, const CNcbiString& ns) {
    return ToBoolean(CNcbiString::Compare(cs, ns) > 0);
}
inline ENcbiBoolean operator>(const CNcbiString& ns1, const CNcbiString& ns2) {
    return ToBoolean(CNcbiString::Compare(ns1, (const Char*)ns2) > 0);
}

inline ENcbiBoolean operator>=(const CNcbiString& ns, const Char *cs) {
    return ToBoolean(CNcbiString::Compare(ns, cs) >= 0);
}
inline ENcbiBoolean operator>=(const Char *cs, const CNcbiString& ns) {
    return ToBoolean(CNcbiString::Compare(cs, ns) >= 0);
}
inline ENcbiBoolean operator>=(const CNcbiString& ns1, const CNcbiString& ns2){
    return ToBoolean(CNcbiString::Compare(ns1, (const Char*)ns2) >= 0);
}


/* CNcbiStringList:: member functions
 */

inline CNcbiString *CNcbiStringList::First() const {
    return (CNcbiString *)CNcbiList::First();
}

inline CNcbiString *CNcbiStringList::Last() const {
    return (CNcbiString *)CNcbiList::Last();
}

inline ENcbiBoolean CNcbiStringList::Remove(const CNcbiString *ns) {
    return CNcbiList::Remove( (void *)ns );
}

inline CNcbiString *CNcbiStringList::RemoveFirst() {
    return (CNcbiString *)CNcbiList::RemoveFirst();
}

inline CNcbiString *CNcbiStringList::RemoveLast() {
    return (CNcbiString *)CNcbiList::RemoveFirst();
}

inline ENcbiBoolean CNcbiStringList::DestroyFirst() {
    CNcbiString *ns = RemoveFirst();
    return ns ? delete ns, kNcbiTrue : kNcbiFalse;
}

inline ENcbiBoolean CNcbiStringList::DestroyLast() {
    CNcbiString *ns = RemoveLast();
    return ns ? delete ns, kNcbiTrue : kNcbiFalse;
}

#endif /* __NCBISTR_HPP__ */
