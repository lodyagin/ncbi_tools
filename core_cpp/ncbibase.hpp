/* $Id: ncbibase.hpp,v 1.2 1997/10/27 18:56:59 vakatov Exp $ */
/*****************************************************************************

 Description: base class for NCBI hierarchy

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

 Classes:
    "CNcbiObject"

 Modification History:
     Sep  2 1997 (grisha) -- originally written

 Bugs and restriction in use:

 Notes:

*****************************************************************************/

#ifndef __NCBIBASE_HPP__
#define __NCBIBASE_HPP__ ncbibase_hpp

#include <ncbidefs.h>
#include <ncbicls.hpp>

/****************************************************************************/
/* Class definitions */
/****************************************************************************/

/* Base class for the whole hierarhy */
class CNcbiObject
{
private:
    static Int4 m_ObjectCount; /* how many objects live in the program */

public:
    /* Constructor/destructor */
    CNcbiObject() { m_ObjectCount++; }
    virtual ~CNcbiObject() { m_ObjectCount--; }

    /* Allocate and return a copy of the object  */
    virtual CNcbiObject* Clone() const = 0;

    static Int4 ObjectCount() { return m_ObjectCount; }
};

#endif /* __NCBIBASE_HPP__ */
