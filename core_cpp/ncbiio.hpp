/* $Id: ncbiio.hpp,v 1.2 1997/10/27 18:57:01 vakatov Exp $ */
/*****************************************************************************

    Description: NCBI file I/O class definition

    Author: Grigoriy Starchenko

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

    Entry Points:

    Modification History:
      Sep 19 1997 - grisha - original written


    Bugs and restriction on use:

    Notes: 

*****************************************************************************/

#ifndef __NCBIIO_HPP__
#define __NCBIIO_HPP__ ncbiio_hpp

/**********************************************************************/
/* Includes */
/**********************************************************************/
#include <ncbibase.hpp>
#include <ncbistr.hpp>

/**********************************************************************/
/* Base class for file I/O */
/**********************************************************************/
class CNcbiIo : public CNcbiObject
{
public:
    /* Constructor/destructor */
    CNcbiIo (void);
    virtual ~CNcbiIo (void);

    /* Member functions */
    const CNcbiString& Name (void) const;
    virtual ENcbiBoolean Open (const CNcbiString&,Int4) = 0;
    virtual ENcbiBoolean Close (void) = 0;
    virtual Int4 Size (void) const;
    virtual Int4 Read (void*,Int4);
    virtual Int4 Write (void*,Int4);
    virtual ENcbiBoolean Gets (CNcbiString*);
    virtual ENcbiBoolean Puts (const CNcbiString&);

protected:
    CNcbiString m_Name;
};

class CNcbiFileMap : public CNcbiIo
{
public:
    /* Constructor/destructor */
    CNcbiFileMap (void);
    virtual ~CNcbiFileMap (void);

    /* Member functions */
    virtual ENcbiBoolean Open (const CNcbiString&,Int4 = 0);
    virtual ENcbiBoolean Close (void);
    virtual Int4 Size (void) const;
    virtual ENcbiBoolean Gets (CNcbiString*);
    void* Image (void) const;
    ENcbiBoolean Begin (void);

private:
    // N/Implemented ==> implemented as a FAKE & hidden!!!
    virtual CNcbiObject* Clone() const { NO_GOOD; return 0; }

private:
    Nlm_MemMapPtr m_Map;
    Int4 m_iCurPos;
};

#endif /* __NCBIIO_HPP__ */
/*EOF*/
