/* ===========================================================================
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
* File Name:  fdlKludge.h
*
* Author:  Jian Ye
*
* Version Creation Date:   10/15/01
*
* $Revision: 6.6 $
*
* File Description:
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: fdlKludge.h,v $
* Revision 6.6  2002/09/11 19:53:09  jianye
* Added url defines
*
* Revision 6.5  2002/08/22 20:32:35  jianye
* add parentheses to bit shift
*
* Revision 6.4  2002/08/21 21:15:32  camacho
* Added #define value for structure link bits
*
* Revision 6.3  2001/10/19 14:40:41  jianye
* *** empty log message ***
*
* Revision 6.2  2001/10/18 19:20:20  jianye
* Initial check in
*
*/

#ifndef FDLKLUDGE_H
#define FDLKLUDGE_H

#include <objloc.h>

#define total_linkout 3
#define linkout_locuslink (1<<0)
#define linkout_unigene   (1<<1)
#define linkout_structure (1<<2)

/* url for linkout*/
#define URL_LocusLink "<a href=\"http://www.ncbi.nlm.nih.gov/LocusLink/list.cgi?Q=%d%s\"><img border=0 height=16 width=16 src=\"/blast/images/L.gif\" alt=\"LocusLink info\"></a>"
#define URL_Unigene "<a href=\"http://www.ncbi.nlm.nih.gov/UniGene/query.cgi?ORG=%s&TEXT=@gi(%d)\"><img border=0 height=16 width=16 src=\"/blast/images/U.gif\" alt=\"UniGene info\"></a>"

#define URL_Structure "<a href=\"http://scarecrow:5701/Structure/cblast_new/cblast.cgi?blast_RID=%s&blast_rep_gi=%d&hit=%d&blast_CD_RID=%s\"><img border=0 height=16 width=16 src=\"/blast/images/R.gif\" alt=\"Structure info\"></a>"

#define URL_Structure_Overview "<a href=\"http://scarecrow:5701/Structure/cblast_new/cblast.cgi?blast_RID=%s&blast_rep_gi=%d&hit=%d&blast_CD_RID=%s\">Structure Linkout Overview</a>"

#endif
