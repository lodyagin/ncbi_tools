/* $Id: ncbiio.cpp,v 1.2 1997/09/22 14:22:33 grisha Exp $ */
/*****************************************************************************

    Description: NCBI file I/O class implementation

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
#ifndef __NCBIIO_CPP__
#define __NCBIIO_CPP__ ncbiio_cpp

/**********************************************************************/
/* Includes */
/**********************************************************************/
#include <ncbiio.hpp>

/**********************************************************************/
/* CNcbiIo implementation */
/**********************************************************************/
/*FCN*/CNcbiIo::CNcbiIo (void)
{
    return;
}                             /* CNcbiIo::CNcbiIo() */

/*FCN*/CNcbiIo::~CNcbiIo (void)
{
    return;
}                             /* CNcbiIo::~CNcbiIo() */

const CNcbiString&
/*FCN*/CNcbiIo::Name (void) const
{
    return m_Name;
}                             /* CNcbiIo::Name() */

Int4
/*FCN*/CNcbiIo::Size (void) const
{
    return -1;
}                             /* CNcbiIo::Size() */

Int4
/*FCN*/CNcbiIo::Read (
  void*,
  Int4
){
    return 0;
}                             /* CNcbiIo::Read() */

Int4
/*FCN*/CNcbiIo::Write (
  void*,
  Int4
){
    return 0;
}                             /* CNcbiIo::Write() */

ENcbiBoolean
/*FCN*/CNcbiIo::Gets (
  CNcbiString*
){
    return kNcbiBad;
}                             /* CNcbiIo::Gets() */

ENcbiBoolean
/*FCN*/CNcbiIo::Puts (
  const CNcbiString&
){
    return kNcbiBad;
}                             /* CNcbiIo::Puts() */

/**********************************************************************/
/* CNcbiFileMap implementation */
/**********************************************************************/
/*FCN*/CNcbiFileMap::CNcbiFileMap (void)
{
    m_Map = NULL;
    m_iCurPos = 0;
}                            /* CNcbiFileMap::CNcbiFileMap() */

/*FCN*/CNcbiFileMap::~CNcbiFileMap (void)
{
    Close();
}                            /* CNcbiFileMap::~CNcbiFileMap() */

ENcbiBoolean
/*FCN*/CNcbiFileMap::Open (
  const CNcbiString& theName,
  Int4
){
    m_Name = theName;
    if ( (m_Map = Nlm_MemMapInit ((CharPtr)m_Name.C_Str())) == NULL ) {
        return kNcbiBad;
    }
    return kNcbiGood;
}                            /* CNcbiFileMap::Open() */

ENcbiBoolean
/*FCN*/CNcbiFileMap::Close (void)
{
    if ( m_Map != NULL ) {
        Nlm_MemMapFini (m_Map);
        m_Map = NULL;
    }
    m_iCurPos = 0;
    return kNcbiGood;
}                            /* CNcbiFileMap::Close() */

Int4
/*FCN*/CNcbiFileMap::Size (void) const
{
    if ( m_Map != NULL ) {
        return m_Map->file_size;
    }
    return 0;
}                            /* CNcbiFileMap::Size() */

void*
/*FCN*/CNcbiFileMap::Image (void) const
{
    if ( m_Map != NULL ) {
        return (void*)m_Map->mmp_begin;
    }
    return NULL;
}                            /* CNcbiFileMap::Image() */

ENcbiBoolean
/*FCN*/CNcbiFileMap::Begin (void)
{
    m_iCurPos = 0;
    return kNcbiGood;
}                            /* CNcbiFileMap::Begin() */

ENcbiBoolean
/*FCN*/CNcbiFileMap::Gets (
  CNcbiString* pOutput
){
    register Int4 iTemp;
    Int4 iNumToSkip;
    register Int4 iMaxSize;
    register CharPtr lpImage;

    if (    (iTemp = m_iCurPos) < (iMaxSize = Size())
         && (lpImage = (CharPtr)Image()) != NULL ) {
        iNumToSkip = 1;             /* by default '\n' */
        while (    iTemp < iMaxSize
                && lpImage[iTemp++] != '\n' ) {
            ;
        }
        if ( iTemp > 1 && lpImage[iTemp-2] == '\r' ) {
            iNumToSkip++;
        }
        pOutput->Clean();
        pOutput->Append (
                     lpImage+m_iCurPos,
                     iTemp-m_iCurPos-iNumToSkip
                 );
        m_iCurPos = iTemp;
    } else {
        return kNcbiBad;
    }

    return kNcbiGood;
}                            /* CNcbiFileMap::Gets() */

#endif /* __NCBIIO_CPP__ */
/*EOF*/
