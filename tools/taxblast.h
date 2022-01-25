/* $Id: taxblast.h,v 6.1 2000/05/17 15:54:39 shavirin Exp $
 * ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* File Name:  $RCSfile: taxblast.h,v $
*
* Authors:  Scott Federhen, Sergei Shavirin
*
* Initial Version Creation Date: 04/04/2000
*
* $Revision: 6.1 $
*
* File Description:
*         Header file for Tax-Blast program
*
* $Log: taxblast.h,v $
* Revision 6.1  2000/05/17 15:54:39  shavirin
* Initial revision in new location.
*
* Revision 1.2  2000/04/26 19:54:03  shavirin
* Added parameter "show_gi" optionaly show gis in seqid in organism list.
*
* Revision 1.1  2000/04/20 14:01:53  shavirin
* Initial revision
*
*
* ==========================================================================
*/

#ifndef TAXBLAST__H
#define TAXBLAST__H

#ifdef __cplusplus
extern "C" {
#endif

    void TXBHtmlReport(SeqAlignPtr sap,   /* output of blast search */
                       FILE *outfile,     /* file where output goes to */ 
                       Boolean query_is_na,  /* If this is DNA sequence? */
                       Boolean db_is_na,     /* If database is DNA? */
                       CharPtr database, /* Database used in BLAST search */
                       CharPtr link_href, /* Link to the regular BLAST results */
                       CharPtr window_name, /* Window name for output */
                       Boolean show_gi); /* Show gis in the organism report */
#ifdef __cplusplus
}
#endif

#endif /*TAXBLAST__H */

#if 0

Int4 CountAligns (SeqAlignPtr sap);
Boolean GetDbMolType (SeqAnnotPtr sap);
HitObjPtr GetAlignData (SeqAlignPtr sap);
Int4 FindTaxid (Int4 taxid, OrgObjPtr orgobj);
OrgObjPtr GetOrgData (HitObjPtr hitobj);
CnamesPtr GetCommonNames (Int4 taxid);
BnamePtr GetBlastName (Int4 taxid, TreePtr tree);
TreePtr GetTreeData (OrgObjPtr orgobj);
LinObjPtr GetLinData (TreePtr tree, Int4 focus);

void TXBPrintReport (FILE *outfile, Uint1 format, Uint1 numhits, Int4 focus,
                     HitObjPtr hitobj, OrgObjPtr orgobj, TreePtr tree);

HitObjPtr HitObjNew (void);
void HitObjFree (HitObjPtr hitobj);
OrgObjPtr OrgObjNew (void);
void OrgObjFree (OrgObjPtr orgobj);
CnamesPtr CnamesNew (void);
void CnamesFree (CnamesPtr cnames);
BnamePtr BnameNew (void);
void BnameFree (BnamePtr bname);
NodeObjPtr NodeObjNew (void);
void NodeObjFree (NodeObjPtr nodeobj);
LinObjPtr LinObjNew (void);
void LinObjFree (LinObjPtr linobj);

void taxreport_maxlen (TreeCursorPtr cursor, Int2 depth, Int2Ptr maxlen);
void taxreport_fn (FILE* outfile, TreeCursorPtr cursor, Int2 depth,
		   Int2 maxlen, Int2 maxhits, Int2 maxorgs);

void traverse_tree (TreePtr tree, void tree_fn());
void traverse_tree_internal (TreeCursorPtr cursor, void tree_fn());
void tree_fnDel (TreeCursorPtr cursor);
ValNodePtr ValNodeDiff (ValNodePtr vnp1, ValNodePtr vnp2);
static int LIBCALLBACK vn_sortfn (VoidPtr ptr1, VoidPtr ptr2);
void ValNodeIntCpy (ValNodePtr PNTR to, ValNodePtr from);
void pad_left (CharPtr line, Int4 depth);
void pad_right (CharPtr line, Int4 max, Char c);
#endif
