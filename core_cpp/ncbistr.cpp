/* $Id: ncbistr.cpp,v 1.11 1997/12/02 21:50:31 vakatov Exp $
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
* Author:  Denis Vakatov
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
*    Test code(use "#define TEST_MODULE_NCBISTR_CPP" to turn it on)
*
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ncbistr.cpp,v $
* Revision 1.11  1997/12/02 21:50:31  vakatov
* CNcbiString::  Added Trim[Right|Left]() and Prune();; now most of aux.
* and editing funcs return reference to the string rather than "void"
*
* Revision 1.10  1997/12/02 00:01:44  vakatov
* Added SetChar(), DeleteChar(), InsertChar() and InsertStr();
* expanded static and member Compare()'s functionality
*
* Revision 1.9  1997/11/12 22:16:49  vakatov
* [CNcbiString, CNcbiBaseStr]  Renamed "Size()" to "BufferSize()"
*
* Revision 1.8  1997/10/28 21:01:15  vakatov
* CNcbiString:: +FreezeBuffer(), +CommitBuffer() -- direct access to the
* string's character buffer
*
* Revision 1.7  1997/10/27 18:57:20  vakatov
* Partially redesigned and fully implemented base NCBI C++
* classes(base, strings, lists, iterators)
*
* ==========================================================================
*/

#include <ncbistr.hpp>


/****************************************************************************
 * CNcbiBaseStr implementation
 ****************************************************************************/

struct SFakeNullStr {
    CNcbiBaseStr::SBigAttr attr;
    Char str[1];
} s_FakeNullStr = {{0, 1, 0, 0}, ""};

CNcbiBaseStr *CNcbiBaseStr::sm_NullStr = (CNcbiBaseStr *)
    &s_FakeNullStr.str[0];

Uint4 CNcbiBaseStr::sm_Counter      = 0;
Uint4 CNcbiBaseStr::sm_BigCounter   = 0;
Uint4 CNcbiBaseStr::sm_SmallCounter = 0;

void *CNcbiBaseStr::operator new(size_t /* not used but must present here */,
                                 const Char *str,
                                 Uint4 str_length,
                                 ESizePolicy sz_policy) {
    Uint4 str_len = 0;
    if ( !str_length ) {
        if ( !str ) {
            sm_Counter--;  // no real new instance created
            return sm_NullStr->Clone();
        }
        str_length = strlen( str );
        str_len = str_length;
    } else if ( str ) {
        str_len = CNcbiString::strnlen(str, str_length);
    }

    Uint4 str_size = str_length + 1;
    if (str_length  &&  sz_policy == eReserveSize)
        str_size = (str_size * 3 / 10 + 2) * 4;

    CNcbiBaseStr *X = NULL;
    if (str_size <= 255) { // small string
        SSmallAttr *A;
        if ( !(A = (SSmallAttr *)new char[sizeof(SSmallAttr) + str_size]) )
            return NULL;
        A->ref_count = 1;
        A->size      = (Uint1)str_size;
        A->length    = (Uint1)str_len;
        X = (CNcbiBaseStr *)(A + 1);
        sm_SmallCounter++;
    } else { // big string
        SBigAttr *A;
        if ( !(A = (SBigAttr *)new char[sizeof(SBigAttr) + str_size]) )
            return NULL;
        A->ref_count = 1;
        A->size      = str_size;
        A->length    = str_len;
        A->is_small  = 0;
        X = (CNcbiBaseStr *)(A + 1);
        sm_BigCounter++;
    }

    *X->Str() = '\0';
    if ( str )
        strncat(X->Str(), str, str_len);
    return X;
}


void CNcbiBaseStr::operator delete(void *ptr) {
    CNcbiBaseStr *X = (CNcbiBaseStr *)ptr;
    if (!X  ||  X->Null())
        return;
    ASSERT ( X->RefCount() );

    if ( X->IsSmall() ) { // small string
        if ( !--X->SmallAttr().ref_count ) {
            ::delete (char *)&X->SmallAttr();
            sm_SmallCounter--;  sm_Counter--;
        }
    } else { // big string
        if (!--X->BigAttr().ref_count  &&  !X->Null()) {
            ::delete (char *)&X->BigAttr();
            sm_BigCounter--;  sm_Counter--;
        }
    }
    ASSERT ( sm_Counter == sm_BigCounter + sm_SmallCounter );
}


CNcbiBaseStr *CNcbiBaseStr::Clone() {
    ASSERT ( Null()  ||  RefCount() );
    if ( IsSmall() ) { // small string
        if (RefCount() == UINT2_MAX)
            return Create( Str() );
        SmallAttr().ref_count++;
    } else { // big string
        ASSERT ( RefCount() != UINT4_MAX );
        BigAttr().ref_count++;
    }
    return this;
}


CNcbiBaseStr *CNcbiBaseStr::Realloc(Uint4       length,
                                    ESizePolicy sz_policy) {
    if (!length  &&  Null())
        return this;  // special case -- Null string

    if (RefCount() == 1  &&  length < BufferSize())
        return this;  // no need to realloc

    CNcbiBaseStr *nbstr = Create(Str(), length, sz_policy);
    if ( !nbstr )
        return NULL;

    Release();
    return nbstr;
}


/****************************************************************************
 * CNcbiString implementation
 ****************************************************************************/
const Char *CNcbiString::kDefaultTrimCharset = " \t\n\r";

void CNcbiString::GetReadyForUpdate(Uint4 length) {
    ASSERT ( m_BS->RefCount() > 0 );
    if (m_BS->RefCount() == 1  &&  length < m_BS->BufferSize())
        return;

    CNcbiBaseStr *prev_BS = m_BS;
    m_BS = CNcbiBaseStr::Create(prev_BS->Str(), length);
    ASSERT ( m_BS );
    prev_BS->Release();
}

CNcbiObject *CNcbiString::Clone() const {
    return new CNcbiString( *this );
}

CNcbiString::~CNcbiString() {
    m_BS->Release();
}

CNcbiString& CNcbiString::DeleteChar(Uint4 pos, Uint4 n_del) {
    Uint4 len = Length();
    if (pos >= len) {
        NO_GOOD;
        return *this;
    }
    if ( !n_del )
        return *this;

    if (pos + n_del >= len) {
        GetReadyForUpdate(pos);
        m_BS->SetLength(pos);
    }
    else {
        len -= n_del;
        GetReadyForUpdate(len);
        Char *s = m_BS->Str() + pos;
        memmove(s, s + n_del, len - pos);
        m_BS->SetLength(len);
    }
    return *this;
}

CNcbiString& CNcbiString::InsertChar(Uint4 pos, Char ch, Uint4 n_ins) {
    Uint4 len = Length();
    if (pos > len) {
        NO_GOOD;
        return *this;
    }
    if ( !n_ins )
        return *this;

    len += n_ins;
    GetReadyForUpdate(len);
    Char *s = m_BS->Str() + pos;
    memmove(s + n_ins, s,  len - pos - n_ins);
    memset (s, ch, n_ins);
    m_BS->SetLength(ch ? len : pos);
    return *this;
}

CNcbiString& CNcbiString::InsertStr(Uint4 pos, const Char *cs, Uint4 n_ins) {
    if (!cs  ||  !*cs  ||  !n_ins)
        return *this;

    Uint4 len = Length();
    if (pos > len) {
        NO_GOOD;
        return *this;
    }
    n_ins = strnlen(cs, n_ins);

    Char *s_old = m_BS->Str();
    CNcbiBaseStr *hold_BS = (s_old <= cs  &&  cs < s_old+BufferSize()) ?
        m_BS->Clone() : 0;  // "cs" is a substring -- hold the allocated buffer
    ASSERT ( !hold_BS  ||  (s_old <= cs  &&  cs <= s_old + len) );

    Uint4 tail_len = len - pos;
    len += n_ins;
    GetReadyForUpdate(len);
    Char *s = m_BS->Str();
    if (hold_BS  &&  s == s_old) {
        s += pos;
        if (s <= cs) {
            memmove(s + n_ins, s,          tail_len);
            memcpy (s,         cs + n_ins, n_ins   );
        }
        else if (s <= cs + n_ins) {
            Uint4 n_ins1 = s - cs;
            Uint4 n_ins2 = n_ins - n_ins1;
            memmove(s + n_ins,  s,          tail_len);
            memcpy (s,          s - n_ins1, n_ins1  );
            memcpy (s + n_ins1, s + n_ins,  n_ins2  );
        }
        else {
            memmove(s + n_ins, s,  tail_len);
            memcpy (s,         cs, n_ins   );
        }
    }
    else {
        s += pos;
        memmove(s + n_ins, s,  tail_len);
        memcpy (s,         cs, n_ins   );
    }
    if ( hold_BS )
        hold_BS->Release();
    m_BS->SetLength(len);
    return *this;
}


CNcbiString& CNcbiString::Append(const CNcbiString& ns) {
    if ( Null() ) {
        return (*this = ns);
    }
    if ( ns.Null() ) {
        return *this;
    }
    if ( Empty() ) {
        return (*this = ns);
    }
    if ( ns.Empty() ) {
        return *this;
    }

    Uint4 len      = ns.Length();
    Uint4 this_len = Length();
    ASSERT ( len  &&  this_len );

    GetReadyForUpdate(this_len + len);
    memcpy(m_BS->Str() + this_len, ns, len);
    m_BS->SetLength(this_len + len);
    return *this;
}

CNcbiString& CNcbiString::Append(const Char *cs, Uint4 len) {
    if (!len || !cs || !*cs)
        return *this;

    len = strnlen(cs, len);
    Uint4 this_len = Length();

    GetReadyForUpdate(this_len + len);
    memcpy(m_BS->Str() + this_len, cs, len);
    m_BS->SetLength(this_len + len);
    return *this;
}

CNcbiString& CNcbiString::Append(Char ch, Uint4 len) {
    if (!ch || !len)
        return *this;

    Uint4 this_len = Length();
    GetReadyForUpdate(this_len + len);
    memset(m_BS->Str() + this_len, ch, len);
    m_BS->SetLength(this_len + len);
    return *this;
}

Uint4 CNcbiString::strnlen(const Char *buf, Uint4 max_len) {
    register Uint4 len = max_len;
    while (*buf  &&  len)  {
        buf++;  len--;
    }
    return (max_len - len);
}

Int4 CNcbiString::Compare(const Char *s1, const Char *s2, Int4 n) {
    if (s1 == s2)
        return  0;
    if ( !s1 )
        return -1;
    if ( !s2 )
        return  1;

    if (n < 0) {
        NO_GOOD;
        n = INT4_MAX; // fool-proof
    }

    const Char *s1_save = s1;
    while (n  &&  *s1  &&  *s1 == *s2) {
        n--; s1++; s2++;
    }
    if (n == 0  ||  *s1 == *s2)
        return 0;

    Int4 pos = s1 - s1_save + 1;
    return (*s1 > *s2) ? pos : -pos;
}


CNcbiString& CNcbiString::ToUpper() {
    GetReadyForUpdate();
    for (Char *str = m_BS->Str();  *str;  str++)
        if ( islower( *str ) )
                *str = (Char)toupper( *str );
    return *this;
}

CNcbiString& CNcbiString::ToLower() {
    GetReadyForUpdate();
    for (Char *str = m_BS->Str();  *str;  str++)
        if ( isupper( *str ) )
            *str = (Char)tolower( *str );
    return *this;
}

static Int4 transliterate_noalnum(const Char *from_str,
                                  Char *to_buf, Int4 len_buf, const Char *fmt)
{
    Int4 lenbuf_1 = len_buf - 1;
    Int4 lenbuf_4 = len_buf - 4;

    register Char ch;
    register Int4 i;

    for (i = 0;  (ch = *from_str) != '\0';  from_str++) {
        if ( isalnum(ch) ) {
            if (i >= lenbuf_1) {
                break;
            }
            to_buf[i++] = ch;
        } else {
            if (i >= lenbuf_4) {
                break;
            }
            i += sprintf(to_buf + i, fmt, (unsigned int)ch);
        }
    }
    to_buf[i] = '\0';

    return i;
}

CNcbiString CNcbiString::UrlVersion() const {
    Char buf[1024*16];
    transliterate_noalnum(C_Str(), buf, sizeof(buf), "%%%02x");
    return buf;
}

CNcbiString CNcbiString::HtmlVersion() const {
    Char buf[1024*16];
    transliterate_noalnum(C_Str(), buf, sizeof(buf), "&#%03d;");
    return buf;
}

CNcbiString& CNcbiString::Trim(const Char *charset) {
    return TrimLeft(charset).TrimRight(charset);
}

CNcbiString& CNcbiString::TrimRight(const Char *charset) {
    if (!charset  ||  !*charset)
        return *this;

    Char *cstr = m_BS->Str();
    if ( !*cstr )
        return *this;

    ASSERT ( Length() > 0 );
    ENcbiBoolean trim_all = kNcbiTrue;
    const Char *str;
    for (str = cstr + Length();  str != cstr; ) {
        str--;
        const Char *s;
        for (s = charset;  *s  &&  *s != *str;  s++);
        if (*s != *str) {
            trim_all = kNcbiFalse;
            break;
        }
    }
    if (str[1] == '\0')  // no trimmable trailing symbols found
        return *this;
    if ( trim_all )      // trim all symbols
        return Clean();

    GetReadyForUpdate();
    m_BS->SetLength(str - cstr + 1);
    return *this;
}

CNcbiString& CNcbiString::TrimLeft(const Char *charset) {
    if (!charset  ||  !*charset)
        return *this;

    Char *cstr = m_BS->Str();
    const Char *str;
    for (str = cstr;  *str;  str++) { // skip trimmable symbols
        const Char *s;
        for (s = charset;  *s  &&  *s != *str;  s++);
        if (*s != *str)
            break;
    }
    if (str == cstr)  // no trimmable leading symbols found
        return *this;
    if ( !*str )      // trim all symbols
        return Clean();

    GetReadyForUpdate();
    Uint4 n_trim = str - cstr;
    ASSERT ( Length() > n_trim );
    Uint4 len = Length() - n_trim;
    memmove(cstr, str, len);
    m_BS->SetLength(len);
    return *this;
}

CNcbiString& CNcbiString::Prune(const Char *charset) {
    if (!charset  ||  !*charset)
        return *this;

    Char *cstr = 0;
    const Char *str;
    for (str = m_BS->Str();  *str;  str++) {
        const Char *s;
        for (s = charset;  *s  &&  *s != *str;  s++);
        if (*s == *str) {
            if ( !cstr ) {
                Uint4 pos = str - m_BS->Str();
                GetReadyForUpdate();
                str = cstr = m_BS->Str() + pos;
            }
        }
        else if ( cstr ) {
            *cstr++ = *str;
        }
    }

    if ( cstr )
        m_BS->SetLength(cstr - m_BS->Str());

    return *this;
}


/****************************************************************************
 * CNcbiStringList implementation
 ****************************************************************************/

CNcbiStringList::~CNcbiStringList() {
    Clean();
}

CNcbiString *CNcbiStringList::Append(const CNcbiString& ns) {
    CNcbiString *_ns = new CNcbiString( ns );
    if ( _ns ) {
        if ( !CNcbiList::Append(_ns) ) {
            delete _ns;
            _ns = NULL;
        }
    }
    return _ns;
}

CNcbiString *CNcbiStringList::Append(const Char* cs) {
    CNcbiString *_ns = new CNcbiString( cs );
    if ( _ns ) {
        if ( !CNcbiList::Append(_ns) ) {
            delete _ns;
            _ns = NULL;
        }
    }
    return _ns;
}

CNcbiString *CNcbiStringList::Prepend(const CNcbiString& ns) {
    CNcbiString *_ns = new CNcbiString( ns );
    if ( _ns ) {
        if ( !CNcbiList::Prepend(_ns) ) {
            delete _ns;
            _ns = NULL;
        }
    }
    return _ns;
}

CNcbiString *CNcbiStringList::Prepend(const Char* cs) {
    CNcbiString *_ns = new CNcbiString( cs );
    if ( _ns ) {
        if ( !CNcbiList::Prepend(_ns) ) {
            delete _ns;
            _ns = NULL;
        }
    }
    return _ns;
}

ENcbiBoolean CNcbiStringList::Destroy(CNcbiString*& ns) {
    if (!ns  ||  !Remove(ns))
        return kNcbiFalse;
    delete ns;
    ns = NULL;
    return kNcbiTrue;
}

void CNcbiStringList::Clean()
{
    while ( DestroyFirst() );
    ASSERT ( IsEmpty() );
}


/****************************************************************************
 * CNcbiString TEST
 ****************************************************************************/
#ifdef TEST_MODULE_NCBISTR_CPP
static void test_str_1()
{ // processing separate strings (all operations but comparisons)
#define ns_arr_SIZE 128
    CNcbiString ns = "012345";
    CNcbiString ns_arr[ns_arr_SIZE];

    Uint4 i;

    for (i = 0;  i < ns_arr_SIZE;  i++) {
        if (i%8 == 0) {
            ns_arr[i] = 'B';
            ns_arr[i] += 'C';
            ns_arr[i] = 'A' + ns_arr[i];
        } else if (i%7 == 0) {
            ns_arr[i] = 'Z' + ns;
            ns_arr[i] += ns_arr[i];
            ns_arr[i] = ns_arr[i] + ns_arr[i] + ns_arr[i];
        } else if (i%6 == 0) {
            ns_arr[i] = ns + "";
        } else if (i%5 == 0) {
            ns_arr[i] = (ns + ns);
        } else if (i%4 == 0) {
            ns_arr[i] = ns;
        } else if (i%3 == 0) {
            ns_arr[i] = ns + '\0';
        } else if (i%2 == 0) {
            ns_arr[i] = (const Char*)0 + ns;
        } else {
            ns_arr[i] = ns + "XXX";
        }
    }

    for (i = 0;  i < ns_arr_SIZE;  i++) {
        if ( !i%5 ) {
            ns_arr[i].Swap( ns );
        } else if ( !i%3 ) {
            ns_arr[i] += 'Y';
        } else if ( !i%7 ) {
            ns_arr[i] += "ZZZ";
        } else {
            ns_arr[i] = 'A' + ns_arr[i];
        }
    }    

    for (i = 1;  i < ns_arr_SIZE;  i++) {
        ns_arr[i-1] = ns_arr[i];
        if ( i%2 )
            ns_arr[i-1].ToLower();
    }

    for (i = 0;  i < ns_arr_SIZE;  i++) {
        if (i%3 == 0)
            ns_arr[i] = (const Char*)NULL;
        else if (i%5 == 0)
            ns_arr[i].Clean();
    }

    CNcbiString nsNull;
    ASSERT ( nsNull.Empty()        &&  nsNull.Null()  &&
             nsNull.Length() == 0  &&  nsNull.BufferSize() == 1 );

    CNcbiString nsEmpty = "";
    ASSERT ( nsEmpty.Empty()        &&  !nsEmpty.Null()  &&
             nsEmpty.Length() == 0  &&  nsEmpty.BufferSize() == 1 );

    CNcbiString nsABC("abc");
    ASSERT ( !nsABC.Empty()       &&  !nsABC.Null()  &&
             nsABC.Length() == 3  &&  nsABC.BufferSize() > 3 );

    ASSERT ( nsNull  + nsNull  == nsNull   );
    ASSERT ( nsEmpty + nsEmpty == nsEmpty  );
    ASSERT ( nsNull  + nsEmpty == nsEmpty  );
    ASSERT ( nsEmpty + nsEmpty == nsEmpty  );
    ASSERT ( nsABC   + nsEmpty == nsABC    );
    ASSERT ( nsEmpty + nsABC   == nsABC    );
    ASSERT ( nsNull  + nsABC   == nsABC    );
    ASSERT ( nsABC   + nsNull  == nsABC    );
    ASSERT ( nsABC   + nsABC   == "abcabc" );

    CNcbiString nsABCABC = nsABC + nsABC;
    ASSERT ( nsABCABC == (nsABC.Swap(nsABC), nsABC += nsABC) );

    CNcbiString ns1 = "xxx";
    CNcbiString ns2 = ns1;
    ns1.SetChar(2, '1');  ASSERT ( ns1 != ns2 );

    ns1 = ns2;  ASSERT ( ns1 == ns2 );
    Char *buf = ns2.FreezeBuffer();
    for (i = 0;  i < ns2.BufferSize()-1;  buf[i++] = 'Z');
    ns2.CommitBuffer(2);  ASSERT ( ns2 == "ZZ" );   ASSERT ( ns1 != ns2 );
    ns2 += 'Y';  ASSERT ( ns2 == "ZZY" );

    ns1.Clean();  ns1 = "1234567890";
    ns2.Clean();  ns2 = "1234567800";
    ASSERT ( ns1.Compare(ns2   ) == 9 );
    ASSERT ( ns1.Compare(ns2,10) == 9 );
    ASSERT ( ns1.Compare(ns2, 8) == 0 );
    Uint4 ui0 = 0;
    ASSERT ( ns1.Compare(ns2, ui0) == 0 );
    ASSERT ( ns1.Compare(ns2) == -ns2.Compare(ns1) );

    ns1.DeleteChar(0   );  ASSERT ( ns1 == "234567890" );
    ns1.DeleteChar(0, 3);  ASSERT ( ns1 == "567890" );
    ns1.DeleteChar(4, 3);  ASSERT ( ns1 == "5678" );
    ns1.DeleteChar(2, 0);  ASSERT ( ns1 == "5678" );
    ns1.DeleteChar(2   );  ASSERT ( ns1 == "568" );

    ns1.InsertChar(1, '7'   );  ASSERT ( ns1 == "5768" );
    ns1.InsertChar(0, '4', 2);  ASSERT ( ns1 == "445768" );
    ns1.InsertChar(6, '9', 3);  ASSERT ( ns1 == "445768999" );

    ns1.InsertStr(0, ns2   );  ASSERT ( ns1 == ns2 + "445768999" );
    ns1.Clean();
    ns1.InsertStr(0, ns2, 8);  ASSERT ( ns1 == "12345678" );
    ns1.InsertStr(0, NULL  );  ASSERT ( ns1 == "12345678" );
    ns1.InsertStr(0, ""    );  ASSERT ( ns1 == "12345678" );
    ns1.InsertStr(0, ns1   );  ASSERT ( ns1 == "1234567812345678" );
    ns1 = "12345";
    ns1.InsertStr(ns1.Length(), ns1, 3);  ASSERT ( ns1 == "12345123" );
    ns1.InsertStr(3, ns1);  ASSERT ( ns1 == "1231234512345123" );

    ns1 = " g\tt\rr\nn" + ns1 + " g\tt\rr\nn";
    ns2 = ns1 + "  \n    ";
    CNcbiString ns3 = "AAA\n" + ns1 + ns2;
    ns1.TrimRight();
    ns1.TrimLeft();        TRACE ( "Trim1:  '%s'\n", (const Char*)ns1 );
    ns2.Trim();            TRACE ( "Trim2:  '%s'\n", (const Char*)ns2 );
    ns1.Prune();           TRACE ( "Prune1: '%s'\n", (const Char*)ns1 );
    ns2.Prune("135\n\r");  TRACE ( "Prune2: '%s'\n", (const Char*)ns2 );
    ns3.Trim("n").Prune("\n\t4").Prune("\r2").Prune("\r2\n\t").Prune("1235B");
    TRACE ( "Prune3: '%s'\n", (const Char*)ns3 );
    ns3 += ns3 + ns3;
    TRACE ( "Add33: '%s'\n", (const Char*)ns3 );
    ns3.Trim().TrimLeft().TrimRight().Prune().Prune("").Prune(0);
    ns3.Clean().Trim().TrimLeft().TrimRight().Prune().Prune("").Prune(0);
}


enum EOper {
    eEQ = 0,
    eNE = 1,
    eLE = 2,
    eGE = 3,
    eLT = 4,
    eGT = 5,

    eN_OPER
};

#define EQ res[eEQ]
#define NE res[eNE]
#define LE res[eLE]
#define GE res[eGE]
#define LT res[eLT]
#define GT res[eGT]

//static const Char *op_str[eN_OPER] = { "==", "!=", "<=", ">=", " <", " >" };


static void test_str_2_res(const ENcbiBoolean res[], int i1, int i2,
                           const Char *cs1, const Char *cs2)
{
    // printout to the logical result table
    TRACE ( "|%2d|%2d|%2d|%2d|%2d|%2d| '%s' '%s'\n",
            (int)EQ, (int)NE, (int)LE, (int)GE, (int)LT, (int)GT,
            cs1 ? cs1 : "<NULL>",  cs2 ? cs2 : "<NULL>" );

    // i1, i2 -dependent rules
    if (i1 == i2)
        ASSERT ( EQ );
    else
        ASSERT ( NE );

    if (i1 < i2)
        ASSERT ( LT );
    if (i1 > i2)
        ASSERT ( GT );

    // generic
    ASSERT ( EQ != NE );
    if ( EQ )
        ASSERT ( GE && LE );
    if ( NE )
        ASSERT ( LT != GT );
    if ( GT )
        ASSERT ( GE && NE );
    if ( LT )
        ASSERT ( LE && NE );
}


static void test_str_2()
{ // processing separate strings (comparison operations)
#define arr_SIZE 6
    static const Char *cs_arr[arr_SIZE] = { // must be in the ascending order
        NULL,
        "",
        "aaa",
        "aaaa",
        "bbb",
        "c"
    };

    CNcbiString ns_arr[arr_SIZE];
    for (int i = 0;  i < arr_SIZE;  i++) {
        ns_arr[i] = cs_arr[i];
    }
    
    for (int i1 = 0;  i1 < arr_SIZE;  i1++) {
        int i2;

        TRACE ( "\n|==|!=|<=|>=| <| >|\n" );

        for (i2 = 0;  i2 < arr_SIZE;  i2++) {
            ENcbiBoolean res[6];
            res[0] = (cs_arr[i1] == ns_arr[i2]);
            res[1] = (cs_arr[i1] != ns_arr[i2]);
            res[2] = (cs_arr[i1] <= ns_arr[i2]);
            res[3] = (cs_arr[i1] >= ns_arr[i2]);
            res[4] = (cs_arr[i1]  < ns_arr[i2]);
            res[5] = (cs_arr[i1]  > ns_arr[i2]);
            test_str_2_res(res, i1, i2, cs_arr[i1], ns_arr[i2]);
        }

        for (i2 = 0;  i2 < arr_SIZE;  i2++) {
            ENcbiBoolean res[6];
            res[0] = (ns_arr[i1] == ns_arr[i2]);
            res[1] = (ns_arr[i1] != ns_arr[i2]);
            res[2] = (ns_arr[i1] <= ns_arr[i2]);
            res[3] = (ns_arr[i1] >= ns_arr[i2]);
            res[4] = (ns_arr[i1]  < ns_arr[i2]);
            res[5] = (ns_arr[i1]  > ns_arr[i2]);
            test_str_2_res(res, i1, i2, ns_arr[i1], ns_arr[i2]);
        }
        for (i2 = 0;  i2 < arr_SIZE;  i2++) {
            ENcbiBoolean res[6];
            res[0] = (ns_arr[i1] == cs_arr[i2]);
            res[1] = (ns_arr[i1] != cs_arr[i2]);
            res[2] = (ns_arr[i1] <= cs_arr[i2]);
            res[3] = (ns_arr[i1] >= cs_arr[i2]);
            res[4] = (ns_arr[i1]  < cs_arr[i2]);
            res[5] = (ns_arr[i1]  > cs_arr[i2]);
            test_str_2_res(res, i1, i2, ns_arr[i1], cs_arr[i2]);
        }
    }
}


static void test_str_3()
{ // engage string lists & iterations
    CNcbiStringList list1;
    for (int i = 1;  i < 128;  i++) {
        if (i%7 == 0) {
            if (i%28 == 0) {
                CNcbiString* pns = list1.Last();
                VERIFY ( list1.Destroy(pns) );  ASSERT ( !pns );
            } else if (i%21 == 0) {
                VERIFY ( list1.DestroyLast() );
            } else if (i%14 == 0) {
                CNcbiString* pns = list1.First();
                VERIFY ( list1.Destroy(pns) );  ASSERT ( !pns );
            } else {
                VERIFY ( list1.DestroyLast() );
            }
        }
        else if (i%3 == 0) {
            Char cs[16];  sprintf(cs, "# %d", i);
            if (i%6 == 0) {
                list1.Append(cs);
            } else {
                CNcbiString *pns = list1.Prepend(cs);  ASSERT ( pns );
                VERIFY ( list1.Destroy(pns) );  ASSERT ( !pns );
            }
        }
        else if (i%2 == 0) {
            if (i%4 == 0) {
                CNcbiString ns("i%4");
                CNcbiString *pns = list1.Append(ns);  ASSERT ( pns );
                VERIFY ( list1.Remove(pns) );
                CNcbiString *pns_dangling = pns;
                pns = list1.Prepend(*pns);  ASSERT ( pns );
                delete pns_dangling;
                VERIFY ( list1.Destroy(pns) );
            } else if ( !list1.IsEmpty() ) {
                CNcbiString ns(*list1.First() + " i%2 ");
                VERIFY ( list1.Prepend(ns + *list1.Last()) );
            }
        }
        else if (i%11 == 0) {
            CNcbiString ns("%else%");
            CNcbiString *pns = list1.Append(ns);
            ASSERT ( pns );
            VERIFY ( list1.Append("{" + *pns + "}") );
            VERIFY ( list1.Destroy( pns ) );
        }
    }


    CNcbiStringListIterator iter1(list1);

    CNcbiStringListIterator iter2(list1);
    iter2.Reset( *iter1.Container() );

    CNcbiStringListIterator* iter3;
    iter3 = (CNcbiStringListIterator*)iter2.Clone();

    TRACE ( "\n\nTotal: %d\n\n", (int)iter3->Container()->Size() );
    for (iter3->Reset();  iter3->IsNotDone();  iter3->Next()) {
        TRACE ( "%s\n", (const char *)(**iter3) );
    }

    delete iter3;
}


static ENcbiBoolean test_str(void (*test_func)())
{ // call the test function;  check for lost string instances
    Uint4 nBaseStr = CNcbiBaseStr::nInstances();

    (*test_func)();

    if (CNcbiBaseStr::nInstances() != nBaseStr) {
        TRACE ( "\n%d %d\n\n\n", CNcbiBaseStr::nInstances(), nBaseStr );
        NO_GOOD;
        return kNcbiBad;
    }
    return kNcbiGood;
}


int main()
{
    if ( !test_str(test_str_1) )
        return 1;

    if ( !test_str(test_str_2) )
        return 2;

    if ( !test_str(test_str_3) )
        return 3;

    return 0;
}
#endif /* TEST_MODULE_NCBISTR_CPP */
