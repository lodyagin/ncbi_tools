/*   cn3dopen.h
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
* File Name:  cn3dopen.h
*
* Author:  Christopher Hogue
*
* Version Creation Date:   1/31/96
*
* $Revision: 6.8 $
*
* File Description: Cn3d file opening routines 
*                   
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* $Log: cn3dopen.h,v $
* Revision 6.8  2000/01/21 15:59:05  lewisg
* add check for binary/ascii files
*
* Revision 6.7  2000/01/04 15:55:51  lewisg
* don't hang on disconnected network and fix memory leak/hang at exit
*
* Revision 6.6  1999/10/29 14:15:30  thiessen
* ran all Cn3D source through GNU Indent to prettify
*
* Revision 6.5  1999/08/04 21:18:01  lewisg
* modularized open operations to allow sequin to launch cn3d
*
* Revision 6.4  1999/01/14 19:07:17  kans
* network availability is configurable
*
* Revision 6.3  1998/06/29 19:28:02  lewisg
* on the fly update of conservation color
*
* Revision 6.2  1998/04/28 19:38:41  lewisg
* codewarrior fixes
*
* Revision 6.1  1998/04/28 15:14:31  lewisg
* moved OpenMimeFileWithDeletion to cn3dopen
*
* Revision 6.0  1997/08/25 18:13:40  madden
* Revision changed to 6.0
*
* Revision 5.0  1996/05/28 14:05:44  ostell
* Set to revision 5.0
*
 * Revision 1.1  1996/02/01  18:47:38  kans
 * Initial revision
 *
*
* ==========================================================================
*/

/* cn3dopen.h */

#ifndef _CN3DOPEN2_
#define _CN3DOPEN2_ 1

#ifdef __cplusplus
extern "C" {
#endif
#include <objmime.h>
#define MAX_MDLNO 1000
#define PRINT_FORM_MIME_NAME "Ncbi-mime-asn1"
#define PRINT_FORM_BIOSTRUC "Biostruc"
extern MenU LIBCALL Cn3D_OpenSub PROTO((MenU m));
extern Boolean OpenMimeFileWithDeletion
    PROTO((CharPtr filename, Boolean removeIt));
extern void LIBCALLBACK fnClearMarkedResidues
    PROTO((PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr));
extern ValNodePtr fnMarkAlignedResidues
    PROTO(
        (PDNMS pdnmsMaster, PDNMS pdnmsSlave,
         BiostrucFeaturePtr pbsfThis));
NLM_EXTERN Boolean MMDB_ReadMime(NcbiMimeAsn1Ptr mime);
NLM_EXTERN void Cn3D_OpenEnd();
NLM_EXTERN void Cn3D_OpenStart();

#ifdef __cplusplus
}
#endif
#endif
