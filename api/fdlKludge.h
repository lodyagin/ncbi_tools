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
* $Revision: 6.18 $
*
* File Description:
*
* Modifications:  
* --------------------------------------------------------------------------
*
*/

#ifndef FDLKLUDGE_H
#define FDLKLUDGE_H

#include <objloc.h>


#define linkout_locuslink (1<<0)
#define linkout_unigene   (1<<1)
#define linkout_structure (1<<2)
#define linkout_geo       (1<<3)
#define linkout_gene      (1<<4)

/* url for linkout*/
#define URL_LocusLink "<a href=\"http://www.ncbi.nlm.nih.gov/LocusLink/list.cgi?Q=%d%s\"><img border=0 height=16 width=16 src=\"/blast/images/L.gif\" alt=\"LocusLink info\"></a>"
#define URL_Unigene "<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=unigene&cmd=search&term=%ld[Nucleotide+UID]\"><img border=0 height=16 width=16 src=\"/blast/images/U.gif\" alt=\"UniGene info\"></a>"

#define URL_Structure "<a href=\"http://www.ncbi.nlm.nih.gov/Structure/cblast/cblast.cgi?blast_RID=%s&blast_rep_gi=%ld&hit=%ld&blast_CD_RID=%s&blast_view=%s&hsp=0&taxname=%s&client=blast\"><img border=0 height=16 width=16 src=\"http://www.ncbi.nlm.nih.gov/Structure/cblast/str_link.gif\" alt=\"Related structures\"></a>"

#define URL_Structure_Overview "<a href=\"http://www.ncbi.nlm.nih.gov/Structure/cblast/cblast.cgi?blast_RID=%s&blast_rep_gi=%ld&hit=%ld&blast_CD_RID=%s&blast_view=%s&hsp=0&taxname=%s&client=blast\">Related Structures</a>"

#define URL_Geo "<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=geo&term=%ld[gi]\"><img border=0 height=16 width=16 src=\"/blast/images/E.gif\" alt=\"Geo\"></a>"
 
#define URL_Gene "<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=search&term=%ld[%s]\"><img border=0 height=16 width=16 src=\"/blast/images/G.gif\" alt=\"Gene info\"></a>"

#endif
