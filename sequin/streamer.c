/*   streamer.c
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
* File Name:  streamer.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   5/17/2005
*
* $Revision: 1.6 $
*
* File Description: 
* This file provides the Main function for the standalone Streaming Editor.
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#ifndef CODECENTER
static char *date_of_compilation = __DATE__;
static char *time_of_compilation = __TIME__;
#else
static char *date_of_compilation = "today";
static char *time_of_compilation = "now";
#endif

#include <vibrant.h>
#include <biosrc.h>
#include <pubdesc.h>
#include <explore.h>
#include <toasn3.h>
#include <tax3api.h>
/* includes for pub lookup */
#include <pmfapi.h>
#include <utilpub.h>
#include <medarch.h>
#include <medutil.h>
#include <mla2api.h>



#ifndef WIN16
CharPtr objPrtMemStr = "PrintTemplateSet ::= {\n" \
"{ name \"StdSeqDesc\" ,\n" \
"format { asn1 \"Seqdesc\" , form block {\n" \
"components {\n" \
"{ asn1 \"Seqdesc.mol-type\" , label \"Molecule type\" , prefix \"\\n\" , form enum { } } ,\n" \
"{ asn1 \"Seqdesc.modif\" , label \"Modifiers\" , prefix \"\\n\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Seqdesc.modif.E\" , form enum { } } } } } ,\n" \
"{ asn1 \"Seqdesc.method\" , label \"Method\" , prefix \"\\n\" , form enum { } } ,\n" \
"{ asn1 \"Seqdesc.name\" , label \"Name\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.title\" , label \"Title\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.org\" , label \"Organism\" , prefix \"\\n\" , form use-template \"StdOrgRef\" } ,\n" \
"{ asn1 \"Seqdesc.comment\" , label \"Comment\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.num\" , label \"Numbering\" , prefix \"\\n\" , form use-template \"StdNumbering\" } ,\n" \
"{ asn1 \"Seqdesc.maploc\" , label \"Map location\" , prefix \"\\n\" , form use-template \"StdDbtag\" } ,\n" \
"{ asn1 \"Seqdesc.pir\" , label \"PIR block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.genbank\" , label \"GenBank block\" , prefix \"\\n\" , form use-template \"StdGBBlock\" } ,\n" \
"{ asn1 \"Seqdesc.pub\" , label \"Citation\" , prefix \"\\n\" , form use-template \"StdPubdesc\" } ,\n" \
"{ asn1 \"Seqdesc.region\" , label \"Region\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.user\" , label \"User Type\" , prefix \"\\n\" , form use-template \"StdUserObj\" } ,\n" \
"{ asn1 \"Seqdesc.sp\" , label \"SWISS-PROT block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.dbxref\" , label \"Cross reference\" , prefix \"\\n\" , form use-template \"StdDbtag\"  } ,\n" \
"{ asn1 \"Seqdesc.embl\" , label \"EMBL block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.create-date\" , label \"Create date\" , prefix \"\\n\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Seqdesc.update-date\" , label \"Update date\" , prefix \"\\n\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Seqdesc.prf\" , label \"PRF block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.pdb\" , label \"PDB block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.het\" , label \"Heterogen\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.source\" , label \"Biological Source\" , prefix \"\\n\" , form use-template \"StdBioSource\" } ,\n" \
"{ asn1 \"Seqdesc.molinfo\" , label \"Molecule Information\" , prefix \"\\n\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"MolInfo.biomol\" , form enum { } } ,\n" \
"{ asn1 \"MolInfo.tech\" , form enum { } } ,\n" \
"{ asn1 \"MolInfo.completeness\" , form enum { } } } } } } } } } ,\n" \
"{ name \"StdSeqFeatLocation\" ,\n" \
"format { asn1 \"Seq-feat.location\" , label \"Location\" , prefix \"\\t\" , form user { printfunc \"StdSeqLocPrint\" } } } ,\n" \
"{ name \"StdSeqFeatProduct\" ,\n" \
"format { asn1 \"Seq-feat.product\" , label \"Product\" , prefix \"\\t\" , form user { printfunc \"StdSeqLocPrint\" } } } ,\n" \
"{ name \"EntrySeqFeatData\" ,\n" \
"labelfrom \"Seq-feat.data\" ,\n" \
"format { asn1 \"Seq-feat.data\" , prefix \"\\t\" , form use-template \"StdSeqFeatData\" } } ,\n" \
"{ name \"StdSeqFeat\" ,\n" \
"labelfrom \"Seq-feat.data\" ,\n" \
"format { asn1 \"Seq-feat\" , prefix \"\\n\" , suffix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.data\" , form use-template \"StdSeqFeatData\" } ,\n" \
"{ asn1 \"Seq-feat\" , form use-template \"StdSeqFeatCommon\" } ,\n" \
"{ asn1 \"Seq-feat.product\" , label \"Product\" , prefix \" \" , form user { printfunc \"StdSeqLocPrint\" } } ,\n" \
"{ asn1 \"Seq-feat.location\" , label \"Location\" , prefix \" \" , form user { printfunc \"StdSeqLocPrint\" } } ,\n" \
"{ asn1 \"Seq-feat.cit\" , label \"Citations\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.cit.pub.E\" , form use-template \"StdPub\" } } } } ,\n" \
"{ asn1 \"Seq-feat.xref\" , label \"Cross-reference\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.xref.E\" , form use-template \"StdSeqFeatXref\" } } } } } } } } ,\n" \
"{ name \"StdSeqFeatData\" ,\n" \
"format { asn1 \"SeqFeatData\" , form block {\n" \
"components {\n" \
"{ asn1 \"SeqFeatData.gene\" , label \"Gene\" , form use-template \"StdGeneRef\" } ,\n" \
"{ asn1 \"SeqFeatData.org\" , label \"Organism\" , form use-template \"StdOrgRef\" } ,\n" \
"{ asn1 \"SeqFeatData.cdregion\" , label \"Coding Region\" , form use-template \"StdCdRegion\" } ,\n" \
"{ asn1 \"SeqFeatData.prot\" , label \"Protein\" , form use-template \"StdProtRef\" } ,\n" \
"{ asn1 \"SeqFeatData.rna\" , label \"RNA\" , form use-template \"StdRNARef\" } ,\n" \
"{ asn1 \"SeqFeatData.pub\" , label \"Citation\" , form use-template \"StdPubdesc\" } ,\n" \
"{ asn1 \"SeqFeatData.seq\" , label \"Sequence\" , form user { printfunc \"StdSeqLocPrint\" } } ,\n" \
"{ asn1 \"SeqFeatData.imp.key\" , label \"Import\" , form use-template \"StdImpFeat\" } ,\n" \
"{ asn1 \"SeqFeatData.region\" , label \"Region\" , form text { } } ,\n" \
"{ asn1 \"SeqFeatData.comment\" , label \"Comment\" , form null NULL } ,\n" \
"{ asn1 \"SeqFeatData.bond\" , label \"Bond\" , form enum { } } ,\n" \
"{ asn1 \"SeqFeatData.site\" , label \"Site\" , form enum { } } ,\n" \
"{ asn1 \"SeqFeatData.rsite\" , label \"Rest. Site\" , form use-template \"StdRsiteRef\" } ,\n" \
"{ asn1 \"SeqFeatData.user\" , label \"User Type\" , form use-template \"StdUserObj\" } ,\n" \
"{ asn1 \"SeqFeatData.txinit\" , label \"TxInit\" , form use-template \"StdTxInit\" } ,\n" \
"{ asn1 \"SeqFeatData.num\" , label \"Numbering\" , form use-template \"StdNumbering\" } ,\n" \
"{ asn1 \"SeqFeatData.psec-str\" , label \"Sec. Struct\" , form enum { } } ,\n" \
"{ asn1 \"SeqFeatData.non-std-residue\" , label \"NonStd Residue\" , form text { } } ,\n" \
"{ asn1 \"SeqFeatData.het\" , label \"Heterogen\" , form text { } } ,\n" \
"{ asn1 \"SeqFeatData.biosrc\" , label \"Biological Source\" , prefix \"\\n\" , form use-template \"StdBioSource\" } } } } } ,\n" \
"{ name \"StdGeneRef\" ,\n" \
"format { asn1 \"Gene-ref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Gene-ref\" , form block {\n" \
"components {\n" \
"{ asn1 \"Gene-ref.locus\" , form text { } } ,\n" \
"{ asn1 \"Gene-ref.allele\" , prefix \" \" , form text { } } } } } ,\n" \
"{ asn1 \"Gene-ref.desc\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Gene-ref.pseudo\" , form boolean {\n" \
"true \"This is a pseudogene.\" } } ,\n" \
"{ asn1 \"Gene-ref.syn\" , label \"Synonyms\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Gene-ref.syn.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Gene-ref.maploc\" , label \"Map Location\" , prefix \" \" , form text { } } ,\n" \
"{ asn1 \"Gene-ref.db\" , label \"Cross Reference\" , prefix \" \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Gene-ref.db.E\" , prefix \"(\" , suffix \")\" , form use-template \"StdDbtag\" } } } } } } } } ,\n" \
"{ name \"StdUserObj\" ,\n" \
"format { asn1 \"User-object\" , label \"User-object\" , prefix \"\\n\" , form block {\n" \
"separator \": \" ,\n" \
"components {\n" \
"{ asn1 \"User-object.class\" , form text { } } ,\n" \
"{ asn1 \"User-object.type\" , form use-template \"StdObjectId\" } } } } } ,\n" \
"{ name \"StdPubOnFeat\" ,\n" \
"format { asn1 \"Pub\" , label \"Citation\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pub\" , form use-template \"StdPub\" } } } } } ,\n" \
"{ name \"StdPub\" ,\n" \
"format { asn1 \"Pub\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pub.gen\" , form use-template \"StdCitGen\" } ,\n" \
"{ asn1 \"Pub.sub\" , form use-template \"StdCitSub\" } ,\n" \
"{ asn1 \"Pub.medline\" , form use-template \"StdMedlineEntry\" } ,\n" \
"{ asn1 \"Pub.muid\" , label \"MEDLINE uid: \" , form text { } } ,\n" \
"{ asn1 \"Pub.pmid\" , label \"PubMed id: \" , form text { } } ,\n" \
"{ asn1 \"Pub.article\" , form use-template \"StdCitArt\" } ,\n" \
"{ asn1 \"Pub.journal\" , form use-template \"StdCitJour\" } ,\n" \
"{ asn1 \"Pub.book\" , form use-template \"StdCitBook\" } ,\n" \
"{ asn1 \"Pub.proc\" , form use-template \"StdCitProc\" } ,\n" \
"{ asn1 \"Pub.patent\" , form use-template \"StdCitPat\" } ,\n" \
"{ asn1 \"Pub.pat-id\" , form use-template \"StdIdPat\" } ,\n" \
"{ asn1 \"Pub.man\" , form use-template \"StdCitLet\" } ,\n" \
"{ asn1 \"Pub.equiv\" , form use-template \"StdPubEquiv\" } } } } } ,\n" \
"{ name \"StdCitGen\" ,\n" \
"format { asn1 \"Cit-gen\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-gen.serial-number\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-gen.date\" , prefix \"(\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Cit-gen.title\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen.cit\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Cit-gen.journal\" , suffix \":\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-gen.issue\" , suffix \";\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen.pages\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdCitSub\" ,\n" \
"format { asn1 \"Cit-sub\" , prefix \"Data Submission \" , form block {\n" \
"components {\n" \
"{ asn1 \"Cit-sub.medium\" , prefix \"on \" , suffix \" \" , form enum { } } ,\n" \
"{ asn1 \"Cit-sub.imp.date\" , prefix \"(\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Cit-sub.authors\" , prefix \"\\n\" , form use-template \"StdAuthList\" } } } } } ,\n" \
"{ name \"StdMedlineEntry\" ,\n" \
"format { asn1 \"Medline-entry\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry\" , form block {\n" \
"separator \"   \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.uid\" , label \"uid\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Medline-entry.em\" , label \"entry month\" , prefix \": \" , form user { printfunc \"StdDatePrint\" } } } } } ,\n" \
"{ asn1 \"Medline-entry.cit\" , form use-template \"StdCitArt\" } ,\n" \
"{ asn1 \"Medline-entry.abstract\" , label \"abstract\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Medline-entry.mesh\" , label \"Mesh Terms\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.mesh.E\" , form block {\n" \
"components {\n" \
"{ asn1 \"Medline-mesh.term\" , form text { } } ,\n" \
"{ asn1 \"Medline-mesh.mp\" , form boolean {\n" \
"true \" (Main Point)\" } } ,\n" \
"{ asn1 \"Medline-mesh.qual\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-mesh.qual.E\" , form block {\n" \
"components {\n" \
"{ asn1 \"Medline-qual.subh\" , form text { } } ,\n" \
"{ asn1 \"Medline-qual.mp\" , form boolean {\n" \
"true \" (Main Point)\" } } } } } } } } } } } } } } ,\n" \
"{ asn1 \"Medline-entry.substance\" , label \"Substance\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.substance.E\" , form block {\n" \
"components {\n" \
"{ asn1 \"Medline-rn.name\" , form text { } } ,\n" \
"{ asn1 \"Medline-rn.type\" , form enum {\n" \
"values {\n" \
"\"\" ,\n" \
"\" CAS: \" ,\n" \
"\"EC \" } } } ,\n" \
"{ asn1 \"Medline-rn.cit\" , form text { } } } } } } } } ,\n" \
"{ asn1 \"Medline-entry.xref\" , label \"Cross Reference\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.xref.E\" , form block {\n" \
"separator \": \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-si.type\" , form enum { } } ,\n" \
"{ asn1 \"Medline-si.cit\" , form text { } } } } } } } } ,\n" \
"{ asn1 \"Medline-entry.gene\" , label \"Possible Gene Symbols\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.gene.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Medline-entry.idnum\" , label \"Support\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.idnum.E\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdCitArt\" ,\n" \
"format { asn1 \"Cit-art\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-art.title\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-art.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-art.from.journal\" , form use-template \"StdCitJour\" } ,\n" \
"{ asn1 \"Cit-art.from.book\" , prefix \"(in) \" , form use-template \"StdCitBook\" } ,\n" \
"{ asn1 \"Cit-art.from.proc\" , prefix \"(in) \" , form use-template \"StdCitProc\" } } } } } ,\n" \
"{ name \"StdCitJour\" ,\n" \
"format { asn1 \"Cit-jour\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Cit-jour.title\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-jour.imp\" , form use-template \"StdImprint\" } } } } } ,\n" \
"{ name \"StdCitBook\" ,\n" \
"format { asn1 \"Cit-book\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-book.title\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-book.coll\" , prefix \"Collection: \" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-book.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-book.imp\" , form use-template \"StdImprint\" } } } } } ,\n" \
"{ name \"StdCitProc\" ,\n" \
"format { asn1 \"Cit-proc\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-proc.book\" , form use-template \"StdCitBook\" } ,\n" \
"{ asn1 \"Cit-proc.meet\" , label \"Meeting \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Meeting.number\" , form text { } } ,\n" \
"{ asn1 \"Meeting.date\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Meeting.place\" , form use-template \"StdAffil\" } } } } } } } } ,\n" \
"{ name \"StdCitPat\" ,\n" \
"format { asn1 \"Cit-pat\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-pat.title\" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-pat\" , form block {\n" \
"components {\n" \
"{ asn1 \"Cit-pat.country\" , suffix \" \" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.doc-type\" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.number\" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.date-issue\" , prefix \" (\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Cit-pat.app-number\" , prefix \" Appl: \" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.app-date\" , prefix \" (\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } } } } } } } } ,\n" \
"{ name \"StdIdPat\" ,\n" \
"format { asn1 \"Id-pat\" , form block {\n" \
"components {\n" \
"{ asn1 \"Id-pat.country\" , suffix \" \" , form text { } } ,\n" \
"{ asn1 \"Id-pat.id.number\" , form text { } } ,\n" \
"{ asn1 \"Id-pat.id.app-number\" , prefix \"Appl: \" , form text { } } } } } } ,\n" \
"{ name \"StdCitLet\" ,\n" \
"format { asn1 \"Cit-let\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-let.type\" , prefix \"[\" , suffix \"]\" , form enum { } } ,\n" \
"{ asn1 \"Cit-let.man-id\" , form text { } } ,\n" \
"{ asn1 \"Cit-let.cit\" , form use-template \"StdCitBook\" } } } } } ,\n" \
"{ name \"StdPubEquiv\" ,\n" \
"format { asn1 \"Pub-equiv\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pub-equiv.E\" , form use-template \"StdPub\" } } } } } ,\n" \
"{ name \"StdTitle\" ,\n" \
"format { asn1 \"Title\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Title.E.trans\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Title.E.name\" , form text { } } ,\n" \
"{ asn1 \"Title.E.tsub\" , form text { } } ,\n" \
"{ asn1 \"Title.E.abr\" , form text { } } ,\n" \
"{ asn1 \"Title.E.iso-jta\" , form text { } } ,\n" \
"{ asn1 \"Title.E.ml-jta\" , label \"MEDLINE\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.jta\" , label \"jta\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.issn\" , label \"ISSN\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.coden\" , label \"CODEN\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.isbn\" , label \"ISBN\" , prefix \": \" , form text { } } } } } } ,\n" \
"{ name \"StdAuthList\" ,\n" \
"format { asn1 \"Auth-list\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Auth-list\" , form user { printfunc \"StdAuthListNamesPrint\" } } ,\n" \
"{ asn1 \"Auth-list.affil\" , form use-template \"StdAffil\" } } } } } ,\n" \
"{ name \"StdAffil\" ,\n" \
"format { asn1 \"Affil\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Affil.str\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.affil\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.div\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.street\" , form text { } } ,\n" \
"{ asn1 \"Affil.std\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Affil.std.city\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.sub\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.country\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdImprint\" ,\n" \
"format { asn1 \"Imprint\" , form block {\n" \
"components {\n" \
"{ asn1 \"Imprint.date\" , prefix \"(\" , suffix \") \" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Imprint.volume\" , form text { } } ,\n" \
"{ asn1 \"Imprint.issue\" , prefix \" (\" , suffix \")\" , form text { } } ,\n" \
"{ asn1 \"Imprint.section\" , prefix \" (\" , suffix \")\" , form text { } } ,\n" \
"{ asn1 \"Imprint.part-sup\" , prefix \" (\" , suffix \")\" , form text { } } ,\n" \
"{ asn1 \"Imprint.pages\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Imprint.prepub\" , prefix \" (\" , suffix \")\" , form enum { } } ,\n" \
"{ asn1 \"Imprint.pub\" , label \"\nPublisher: \" , form use-template \"StdAffil\" } ,\n" \
"{ asn1 \"Imprint.cprt\" , label \" Copyright: \" , form user { printfunc \"StdDatePrint\" } } } } } } ,\n" \
"{ name \"StdSeqFeatXref\" ,\n" \
"format { asn1 \"SeqFeatXref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"SeqFeatXref.id\" , label \"Id=\" , form use-template \"StdFeatId\" } ,\n" \
"{ asn1 \"SeqFeatXref.data\" , form use-template \"StdSeqFeatData\" } } } } } ,\n" \
"{ name \"StdOrgRef\" ,\n" \
"format { asn1 \"Org-ref\" , label \"Org-ref\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.taxname\" , form text { } } ,\n" \
"{ asn1 \"Org-ref.common\" , prefix \"(\" , suffix \")\" , form text { } } } } } ,\n" \
"{ asn1 \"Org-ref.mod\" , label \"Modifiers\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.mod.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Org-ref.db\" , label \"Cross Reference\" , prefix \" \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.db.E\" , prefix \"(\" , suffix \")\" , form use-template \"StdDbtag\" } } } } ,\n" \
"{ asn1 \"Org-ref.syn\" , label \"Synonyms\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.syn.E\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdBioSource\" ,\n" \
"format { asn1 \"BioSource\" , label \"BioSource\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"BioSource.genome\" , form enum { } } ,\n" \
"{ asn1 \"BioSource.org\" , label \"Organism\" , form use-template \"StdOrgRef\" } } } } } ,\n" \
"{ name \"StdCdRegion\" ,\n" \
"format { asn1 \"Cdregion\" , label \"Cdregion\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cdregion.orf\" , form boolean {\n" \
"true \"Uncharacterized Open Reading Frame\" } } ,\n" \
"{ asn1 \"Cdregion.frame\" , label \"Reading Frame = \" , form enum { } } ,\n" \
"{ asn1 \"Cdregion.code\" , label \"Genetic Code: \" , suffix \";\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Genetic-code.E.name\" , form text { } } ,\n" \
"{ asn1 \"Genetic-code.E.id\" , label \"id= \" , form text { } } } } } ,\n" \
"{ asn1 \"Cdregion.conflict\" , form boolean {\n" \
"true \"Translation conflicts with protein sequence\" } } ,\n" \
"{ asn1 \"Cdregion.stops\" , prefix \"Translation contains \" , suffix \" stop codons\" , form text { } } ,\n" \
"{ asn1 \"Cdregion.gaps\" , prefix \"Translation contains \" , suffix \" gaps when aligned to protein\" , form text { } } ,\n" \
"{ asn1 \"Cdregion.mismatch\" , prefix \"Translation contains \" , suffix \" mismatches when aligned to protein\" , form text { } } } } } } ,\n" \
"{ name \"StdProtRef\" ,\n" \
"format { asn1 \"Prot-ref\" , label \"Prot-ref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.name\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.name.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Prot-ref.desc\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Prot-ref.processed\" , form enum { } } ,\n" \
"{ asn1 \"Prot-ref.ec\" , label \"ec\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.ec.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Prot-ref.activity\" , label \"activity\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.activity.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Prot-ref.db\" , label \"Cross Reference\" , prefix \" \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.db.E\" , prefix \"(\" , suffix \")\" , form use-template \"StdDbtag\" } } } } } } } } ,\n" \
"{ name \"StdRNARef\" ,\n" \
"format { asn1 \"RNA-ref\" , label \"RNA-ref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"RNA-ref.type\" , form enum { } } ,\n" \
"{ asn1 \"RNA-ref.pseudo\" , form boolean {\n" \
"true \"This is an RNA pseudogene.\" } } ,\n" \
"{ asn1 \"RNA-ref.ext.name\" , form text { } } } } } } ,\n" \
"{ name \"StdPubdesc\" ,\n" \
"format { asn1 \"Pubdesc\" , label \"Pubdesc\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pubdesc.pub\" , form use-template \"StdPubEquiv\" } ,\n" \
"{ asn1 \"Pubdesc\" , prefix \"In this article:\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pubdesc.name\" , label \"name=\" , form text { } } ,\n" \
"{ asn1 \"Pubdesc.fig\" , label \"figure=\" , form text { } } ,\n" \
"{ asn1 \"Pubdesc.poly-a\" , form boolean {\n" \
"true \"poly(A) shown\" } } ,\n" \
"{ asn1 \"Pubdesc.maploc\" , label \"map location=\" , form text { } } ,\n" \
"{ asn1 \"Pubdesc.num\" , form use-template \"StdNumbering\" } ,\n" \
"{ asn1 \"Pubdesc.numexc\" , form boolean {\n" \
"true \"numbering inconsistent\" } } } } } ,\n" \
"{ asn1 \"Pubdesc.comment\" , form text { } } } } } } ,\n" \
"{ name \"StdImpFeat\" ,\n" \
"format { asn1 \"Imp-feat.key\" , label \"Imp-feat\" , form text { } } } ,\n" \
"{ name \"StdRsiteRef\" ,\n" \
"format { asn1 \"Rsite-ref\" , label \"Rsite-ref\" , form block {\n" \
"components {\n" \
"{ asn1 \"Rsite-ref.str\" , form text { } } ,\n" \
"{ asn1 \"Rsite-ref.std\" , form use-template \"StdDbtag\" } } } } } ,\n" \
"{ name \"StdTxInit\" ,\n" \
"format { asn1 \"Txinit\" , label \"TxInit\" , form block {\n" \
"components {\n" \
"{ asn1 \"Txinit.name\" , form text { } } } } } } ,\n" \
"{ name \"StdNumbering\" ,\n" \
"format { asn1 \"Numbering\" , label \"Numbering\" , form null NULL } } ,\n" \
"{ name \"StdGBBlock\" ,\n" \
"format { asn1 \"GB-block\" , label \"GenBank-block\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"GB-block.extra-accessions\" , label \"Extra accessions\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"GB-block.extra-accessions.E\" , form text { } } } } } ,\n" \
"{ asn1 \"GB-block.keywords\" , label \"Keywords\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"GB-block.keywords.E\" , form text { } } } } } ,\n" \
"{ asn1 \"GB-block.source\" , label \"Source: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.origin\" , label \"Origin: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.div\" , label \"Division: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.taxonomy\" , label \"Taxonomy: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.date\" , label \"Date: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.entry-date\" , label \"Entry date: \" , form user { printfunc \"StdDatePrint\" } } } } } } ,\n" \
"{ name \"StdFeatId\" ,\n" \
"format { asn1 \"Feat-id\" , form block {\n" \
"components {\n" \
"{ asn1 \"Feat-id.gibb\" , label \"GenInfo Backbone: \" , form text { } } ,\n" \
"{ asn1 \"Feat-id.giim.id\" , label \"GenInfo Import Id: \" , form text { } } ,\n" \
"{ asn1 \"Feat-id.local\" , label \"Local: \" , form use-template \"StdObjectId\" } ,\n" \
"{ asn1 \"Feat-id.general\" , form use-template \"StdDbtag\" } } } } } ,\n" \
"{ name \"StdSeqFeatCommon\" ,\n" \
"format { asn1 \"Seq-feat\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.id\" , label \"Id=\" , form use-template \"StdFeatId\" } ,\n" \
"{ asn1 \"Seq-feat.title\" , form text { } } ,\n" \
"{ asn1 \"Seq-feat\" , suffix \";\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.partial\" , form boolean {\n" \
"true \"Partial\" } } ,\n" \
"{ asn1 \"Seq-feat.except\" , form boolean {\n" \
"true \"Biological Exception\" } } ,\n" \
"{ asn1 \"Seq-feat.exp-ev\" , label \"Evidence\" , prefix \" is \" , form enum { } } } } } ,\n" \
"{ asn1 \"Seq-feat.comment\" , form text { } } ,\n" \
"{ asn1 \"Seq-feat.ext\" , form use-template \"StdUserObj\" } ,\n" \
"{ asn1 \"Seq-feat.qual\" , label \"Qualifiers\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.qual.E\" , prefix \"/\" , form block {\n" \
"separator \"= \" ,\n" \
"components {\n" \
"{ asn1 \"Gb-qual.qual\" , form text { } } ,\n" \
"{ asn1 \"Gb-qual.val\" , form text { } } } } } } } } } } } } ,\n" \
"{ name \"StdDbtag\" ,\n" \
"format { asn1 \"Dbtag\" , form block {\n" \
"components {\n" \
"{ asn1 \"Dbtag.db\" , suffix \": \" , form text { } } ,\n" \
"{ asn1 \"Dbtag.tag\" , form use-template \"StdObjectId\" } } } } } ,\n" \
"{ name \"StdObjectId\" ,\n" \
"format { asn1 \"Object-id\" , form block {\n" \
"components {\n" \
"{ asn1 \"Object-id.id\" , form text { } } ,\n" \
"{ asn1 \"Object-id.str\" , form text { } } } } } } };\n";
#else
CharPtr objPrtMemStr = "";
#endif

Boolean  useTaxon = FALSE;

static PubdescEditProcs    pubedprocs;
static BioSourceEditProcs  biosrcedprocs;
Boolean  useEntrez = FALSE;
Boolean  indexerVersion = FALSE;


static Int2  taxonCount;

static Int4 TaxLookup (SeqEntryPtr sep, Boolean strip, Boolean correct,
                              MonitorPtr mon)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Int4          rsult;
  Char          str [32];

  rsult = 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        rsult += TaxLookup (sep, strip, correct, mon);
      }
      return rsult;
    }
  }
  if (mon != NULL) {
    taxonCount++;
    sprintf (str, "Processing Component %d", (int) taxonCount);
    MonitorStrValue (mon, str);
  }

  oldscope = SeqEntrySetScope (sep);
  rsult = SeqEntryToAsn3Ex (sep, strip, correct, TRUE, NULL, Tax3MergeSourceDescr, FALSE, FALSE);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, sep);
  SeqEntrySetScope (oldscope);
  return rsult;
}

static Boolean LookupTaxonomyFunc (Uint2 entityID)

{
  SeqEntryPtr  sep;
  MonitorPtr   mon;
  Int4        rsult;
  ErrSev      sev;

  if (entityID < 1) return FALSE;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return FALSE;
  
  rsult = 0;
  sev = ErrSetMessageLevel (SEV_FATAL);

  taxonCount = 0;
  WatchCursor ();
  mon = MonitorStrNewEx ("Taxonomy Lookup", 40, FALSE);
  MonitorStrValue (mon, "Processing Organism Info");
  Update ();

  ObjMgrSetDirtyFlag (entityID, TRUE);

  EntryMergeDupBioSources (sep); /* do before and after SE2A3 */

  Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
  rsult = TaxLookup (sep, TRUE, FALSE, mon);
  MonitorStrValue (mon, "Closing Taxon");
  Update ();
  MonitorFree (mon);
  ArrowCursor ();
  Update ();

  ErrSetMessageLevel (sev);
  ErrClear ();
  ErrShow ();

  return TRUE;
}

static Int2 LclGetSequinAppParam (CharPtr section, CharPtr type, CharPtr dflt, CharPtr buf, Int2 buflen)

{
  Int2  rsult;

  rsult = GetAppParam ("SEQUINCUSTOM", section, type, NULL, buf, buflen);
  if (rsult) return rsult;
  rsult = GetAppParam ("SEQUIN", section, type, dflt, buf, buflen);
  return rsult;
}

/* section for enabling pubmed lookups */
static Boolean log_mla_asn = FALSE;
static Boolean log_mla_set = FALSE;

static ValNodePtr LookupAnArticleFuncNew (ValNodePtr oldpep, BoolPtr success)

{
  CitArtPtr      cap = NULL;
  ArticleIdPtr   ids;
  MlaBackPtr     mbp;
  MlaRequestPtr  mrp;
  Int4           pmid = 0;
  ValNodePtr     pub = NULL;
  ValNodePtr     vnp;
#ifdef OS_UNIX
  CharPtr        str;
#endif

  if (success != NULL) {
    *success = FALSE;
  }
  if (oldpep == NULL) return NULL;

#ifdef OS_UNIX
  if (! log_mla_set) {
    str = (CharPtr) getenv ("LOG_MLA_ASN");
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "TRUE") == 0) {
        log_mla_asn = TRUE;
      }
    }
    log_mla_set = TRUE;
  }
#endif

  for (vnp = oldpep; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == PUB_PMid) {
      pmid = (Int4) vnp->data.intvalue;
    }
  }

  if (cap != NULL) {
    mrp = Mla2CreateCitArtMatchRequest (cap);
    if (mrp != NULL) {
      if (log_mla_asn) {
        LaunchAsnTextViewer ((Pointer) mrp, (AsnWriteFunc) MlaRequestAsnWrite, "citart match request");
      }
      mbp = Mla2SynchronousQuery (mrp);
      mrp = MlaRequestFree (mrp);
      if (mbp != NULL) {
        if (log_mla_asn) {
          LaunchAsnTextViewer ((Pointer) mbp, (AsnWriteFunc) MlaBackAsnWrite, "citart match result");
        }
        pmid = Mla2ExtractCitMatchReply (mbp);
        mbp = MlaBackFree (mbp);
      }
    }
  }

  if (pmid > 0) {
    mrp = Mla2CreatePubFetchRequest (pmid);
    if (mrp != NULL) {
      if (log_mla_asn) {
        LaunchAsnTextViewer ((Pointer) mrp, (AsnWriteFunc) MlaRequestAsnWrite, "get pub request");
      }
      mbp = Mla2SynchronousQuery (mrp);
      mrp = MlaRequestFree (mrp);
      if (mbp != NULL) {
        if (log_mla_asn) {
          LaunchAsnTextViewer ((Pointer) mbp, (AsnWriteFunc) MlaBackAsnWrite, "get pub result");
        }
        cap = Mla2ExtractPubFetchReply (mbp);
        if (cap != NULL) {
          if (log_mla_asn) {
            LaunchAsnTextViewer ((Pointer) cap, (AsnWriteFunc) CitArtAsnWrite, "before");
          }
          ChangeCitArtMLAuthorsToSTD (cap);
          if (log_mla_asn) {
            LaunchAsnTextViewer ((Pointer) cap, (AsnWriteFunc) CitArtAsnWrite, "after");
          }
          for (ids = cap->ids; ids != NULL; ids = ids->next) {
            if (ids->choice != ARTICLEID_PUBMED) continue;
            if (ids->data.intvalue != pmid) {
              Message (MSG_POSTERR, "CitArt ID %ld does not match PMID %ld",
                       (long) ids->data.intvalue, (long) pmid);
            }
          }
          pub = ValNodeAddPointer (NULL, PUB_Article, (Pointer) cap);
          ValNodeAddInt (&pub, PUB_PMid, pmid);
          if (success != NULL) {
            *success = TRUE;
          }
        }
        mbp = MlaBackFree (mbp);
      }
    }
  }

  return pub;
}


static Boolean HasMoreThanOneISOJTA (TitleMsgPtr titles)
{
  TitleMsgPtr tmp;
  ValNodePtr  ttl;
  Int4 num_found = 0;

  for (tmp = titles; tmp != NULL; tmp = tmp->next) {
    for (ttl = tmp->title; ttl != NULL; ttl = ttl->next) {
      if (ttl->choice == Cit_title_iso_jta) {
        if (num_found > 0) {
          return TRUE;
        }
        num_found++;
        break;
      }
    }
  }
  return FALSE;
}


static CharPtr GetExtendedDescription (TitleMsgPtr tmp) 
{
  CharPtr iso_jta = NULL, name = NULL, issn = NULL, extended = NULL;
  ValNodePtr ttl;

  if (tmp == NULL) {
    return NULL;
  }
  
  for (ttl = tmp->title; ttl != NULL; ttl = ttl->next) {
    if (ttl->choice == Cit_title_iso_jta) {
      iso_jta = ttl->data.ptrvalue;
    } else if (ttl->choice == Cit_title_name) {
      name = ttl->data.ptrvalue;
    } else if (ttl->choice == Cit_title_issn) {
      issn = ttl->data.ptrvalue;
    }
  }
  if (iso_jta == NULL) {
    return NULL;
  }
  if (name == NULL && issn == NULL) {
    extended = StringSave (iso_jta);
  } else if (name == NULL) {
    extended = (CharPtr) MemNew (sizeof (Char) * (StringLen (iso_jta) + StringLen (issn) + 5));
    sprintf (extended, "%s||(%s)", iso_jta, issn);
  } else if (issn == NULL) {
    extended = (CharPtr) MemNew (sizeof (Char) * (StringLen (iso_jta) + StringLen (name) + 5));
    sprintf (extended, "%s||(%s)", iso_jta, name);
  } else {
    extended = (CharPtr) MemNew (sizeof (Char) * (StringLen (iso_jta) + StringLen (name) + StringLen (issn) + 6));
    sprintf (extended, "%s||(%s:%s)", iso_jta, name, issn);
  }
  return extended;
}


static Boolean LookupJournalFuncNew (CharPtr title, size_t maxsize, Int1Ptr jtaType, ValNodePtr PNTR all_titlesP)

{
  ValNodePtr       first = NULL;
  ValNodePtr       last = NULL;
  MlaBackPtr       mbp;
  MlaRequestPtr    mrp;
  CharPtr          str, end;
  TitleMsgListPtr  tlp;
  TitleMsgPtr      tmp;
  ValNodePtr       ttl;
  ValNodePtr       vnp;

  if (all_titlesP != NULL) {
    *all_titlesP = NULL;
  }
  if (jtaType != NULL) {
    *jtaType = 0;
  }
  if (StringHasNoText (title)) return FALSE;

#ifdef OS_UNIX
  if (! log_mla_set) {
    str = (CharPtr) getenv ("LOG_MLA_ASN");
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "TRUE") == 0) {
        log_mla_asn = TRUE;
      }
    }
    log_mla_set = TRUE;
  }
#endif

  WatchCursor ();
  Update ();

  mrp = Mla2CreateJournalTitleRequest (title);
  if (mrp != NULL) {
  if (log_mla_asn) {
    LaunchAsnTextViewer ((Pointer) mrp, (AsnWriteFunc) MlaRequestAsnWrite, "lookup journal request");
  }
    mbp = Mla2SynchronousQuery (mrp);
    mrp = MlaRequestFree (mrp);
    if (mbp != NULL) {
      if (log_mla_asn) {
        LaunchAsnTextViewer ((Pointer) mbp, (AsnWriteFunc) MlaBackAsnWrite, "lookup journal result");
      }
      tlp = Mla2ExtractJournalTitleReply (mbp);
      if (tlp != NULL) {
        if (HasMoreThanOneISOJTA (tlp->titles)) {
          for (tmp = tlp->titles; tmp != NULL; tmp = tmp->next) {
            str = GetExtendedDescription (tmp);
            if (str != NULL) {
              ValNodeAddPointer (&first, 0, str);
            }
          }
        } else {
          for (tmp = tlp->titles; tmp != NULL; tmp = tmp->next) {
            for (ttl = tmp->title; ttl != NULL; ttl = ttl->next) {
              if (ttl->choice != Cit_title_iso_jta) continue;
              str = (CharPtr) ttl->data.ptrvalue;
              if (StringHasNoText (str)) continue;
              vnp = ValNodeCopyStr (&last, ttl->choice, str);
              if (first == NULL) {
                first = vnp;
              }
              last = vnp;
            }
          }
        }
        if (first != NULL) {
          str = (CharPtr) first->data.ptrvalue;
          StringNCpy_0 (title, str, maxsize);
          end = StringSearch (title, "||");
          if (end != NULL) {
            *end = 0;
          }
          if (jtaType != NULL) {
            *jtaType = Cit_title_iso_jta;
          }
        }
        if (all_titlesP == NULL) {
          first = ValNodeFreeData (first);
        } else {
          *all_titlesP = first;
        }
        tlp = TitleMsgListFree (tlp);
      }
      mbp = MlaBackFree (mbp);
    }
  }

  ArrowCursor ();
  Update ();

  return TRUE;
}


static Boolean LoadReaderData ()
{
  Char  str [64];

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return FALSE;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return FALSE;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return FALSE;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return FALSE;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return FALSE;
  }
  
  if (! PrintTemplateSetLoadEx ("objprt.prt", objPrtMemStr)) {
    ArrowCursor ();
    Message (MSG_FATAL, "PrintTemplateSetLoad objprt.prt failed");
    return 0;
  }

  
  SetupGeneticCodes ();

  if (! LoadOrganismTable ()) {
    ArrowCursor ();
    Message (MSG_POSTERR, "LoadOrganismTable failed");
  }
  
  MemSet ((Pointer) (&pubedprocs), 0, sizeof (PubdescEditProcs));
  pubedprocs.lookupArticle = LookupAnArticleFuncNew;
  pubedprocs.lookupJournal = LookupJournalFuncNew;
  SetAppProperty ("PubdescEditForm", &pubedprocs);

  useTaxon = FALSE;
 #ifdef USE_TAXON
  useTaxon = TRUE;
#endif
  if (LclGetSequinAppParam ("SETTINGS", "USETAXON", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useTaxon = TRUE;
    }
  }
  
  if (LclGetSequinAppParam ("SETTINGS", "INDEXERVERSION", NULL, str, sizeof (str))) {
      SetAppProperty ("InternalNcbiSequin", (void *) 1024);
  }
#ifdef INTERNAL_NCBI_SEQUIN
  SetAppProperty ("InternalNcbiSequin", (void *) 1024);
#endif

 
  MemSet ((Pointer) (&biosrcedprocs), 0, sizeof (BioSourceEditProcs));
  if (useTaxon) {
    biosrcedprocs.lookupTaxonomy = LookupTaxonomyFunc;
  }
  SetAppProperty ("BioSourcEditForm", &biosrcedprocs);

  return TRUE;
}

static void CloseStreamEditorForm (WindoW w)
{
  if (Message (MSG_OKC, "Are you sure you want to quit?") == ANS_CANCEL) {
    return;
  }
  Remove (w);
  QuitProgram ();
}


typedef struct streamfilterform {
  FORM_MESSAGE_BLOCK

  DialoG pub_dlg;

  CharPtr    input_file;
  CharPtr    output_file;
  ValNodePtr desc_stream_list;
  SeqIdPtr   sip_list;
  Boolean    is_binary;
  Boolean    is_batch;
  Boolean    is_submit;
} StreamFilterFormData, PNTR StreamFilterFormPtr;


static void AcceptFilteringChanges (ButtoN b)
{
  StreamFilterFormPtr frm;
  FILE                *fp;
  FILE                *outfp;
  Char                buf[PATH_MAX];
  Boolean             move_output = FALSE;

  frm = (StreamFilterFormPtr) GetObjectExtra (b);
  if (frm == NULL) 
  {
    return;
  }

  fp = FileOpen (frm->input_file, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", frm->input_file);
    return;
  }
  if (StringHasNoText (frm->output_file) || StringCmp (frm->input_file, frm->output_file) == 0) {
    sprintf (buf, "%s.filter_tmp", frm->input_file);
    outfp = FileOpen (buf, "w");       
    move_output = TRUE;
  } else {
    outfp = FileOpen (frm->output_file, "w");
  }
  WriteAsnWithReplacedDescriptors (frm->desc_stream_list, fp, outfp, frm->is_binary, frm->is_batch, frm->is_submit);
  FileClose (fp);
  FileClose (outfp);  

  if (move_output) {
    FileRemove (frm->input_file);
    FileRename (buf, frm->input_file);
  }

  Remove ((WindoW) frm->form);
  QuitProgram ();
}


static void CancelStreamEditor (ButtoN b)
{
  StreamFilterFormPtr frm;

  frm = (StreamFilterFormPtr) GetObjectExtra (b);
  if (frm == NULL) 
  {
    return;
  }

  if (Message (MSG_OKC, "Are you sure you want to quit?") == ANS_CANCEL) {
    return;
  }

  Remove ((WindoW) frm->form);
  QuitProgram ();
}


static void CleanupStreamFilteringForm (GraphiC g, VoidPtr data)

{
  StreamFilterFormPtr  frm;

  frm = (StreamFilterFormPtr) data;
  if (frm != NULL) {
    frm->desc_stream_list = DescStreamListFree(frm->desc_stream_list);
  }
  MemFree (data);
}


static WindoW MakeStreamFilteringWindow (CharPtr input_file, CharPtr output_file, Boolean is_binary, Boolean is_batch, Boolean is_submit)
{
  FILE         *fp;
  ValNodePtr   tmp;
  WindoW       w;
  GrouP        h, g, c;
  ButtoN       b;
  StreamFilterFormPtr frm;
  SeqIdPtr     sip_list = NULL;
  time_t       t1, t2;

  fp = FileOpen (input_file, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", input_file);
    return NULL;
  }

  t1 = time (NULL);
  tmp = StreamAsnForDescriptors(fp, is_binary, is_batch, is_submit, &sip_list);
  t2 = time (NULL);
  FileClose (fp);

  frm = (StreamFilterFormPtr) MemNew (sizeof (StreamFilterFormData));

  w = FixedWindow (-50, -33, -10, -10, "Pub Editing", CloseStreamEditorForm);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  frm->form = (ForM) w;
  SetObjectExtra (w, frm, CleanupStreamFilteringForm);
  
  frm->input_file = input_file;
  frm->output_file = output_file;
  frm->desc_stream_list = tmp;
  frm->sip_list = sip_list;
  frm->is_binary = is_binary;
  frm->is_batch = is_batch;
  frm->is_submit = is_submit;

  g = HiddenGroup (h, -1, 0, NULL);
  frm->pub_dlg = DescriptorStreamEditor (g, NULL, NULL);
  SetDescriptorStreamEditorIdList(frm->pub_dlg, frm->sip_list);
  PointerToDialog (frm->pub_dlg, frm->desc_stream_list);
  t1 = time (NULL);

  c = HiddenGroup (h, 2, 0, NULL);
  b = DefaultButton (c, "Accept", AcceptFilteringChanges);
  SetObjectExtra (b, frm, NULL);
  b = PushButton (c, "Cancel", CancelStreamEditor);
  SetObjectExtra (b, frm, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  return w;
}


#define STREAMER_APP_VER "1.0"

CharPtr STREAMER_APPLICATION = STREAMER_APP_VER;


Int2 Main (void)

{

  WindoW       w;
  Char         app [64];
  CharPtr      input_file = NULL;
  CharPtr      output_file = NULL;
  Boolean      is_binary = FALSE, is_batch = FALSE, is_submit = FALSE;
  Int2         i;
  Int4         argc = GetArgc();
  CharPtr PNTR argv = GetArgv();
  CharPtr      type = NULL;

#if defined(OS_MAC) && !defined(OS_UNIX_DARWIN)
  long           sysVer;
#endif

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ProcessUpdatesFirst (FALSE);
  
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

#if defined(OS_MAC) && !defined(OS_UNIX_DARWIN)
  if ( Gestalt (gestaltSystemVersion, &sysVer) == noErr) {
    /* system version in low order word is hexadecimal */
    if (sysVer >= 4096) {
      Message (MSG_OK, "You are running on MacOS X and should use the native version of Sequin, not SequinOS9");
    }
  }
#endif

  if (! LoadReaderData ())
  {
    QuitProgram ();
  }  

  sprintf (app, "streamer %s", STREAMER_APPLICATION);
  for (i = 1;  i < argc;  i++)
  {
    if (StringCmp (argv[i], "-in") == 0)
    {
      input_file = argv[i + 1];
      i++;
    }
    else if (StringCmp (argv[i], "-out") == 0)
    {
      output_file = argv[i + 1];
      i++;
    }
    else if (StringCmp (argv[i], "-b") == 0) 
    {
      is_binary = TRUE;
    }
    else if (StringCmp (argv[i], "-a") == 0)
    {
      if (i + 1 < argc) 
      {
        type = argv[i + 1];
        i++;
        if (StringCmp (type, "n") == 0) 
        {
          is_batch = TRUE;
        }
        else if (StringCmp (type, "m") == 0)
        {
          is_submit = TRUE;
        }
      }
    }
    else
    {
      Message (MSG_ERROR, "Unrecognized argument for %s", app);
      QuitProgram();
    }
  }

  if (StringHasNoText (input_file)) {
    Message (MSG_ERROR, "No input file specified.");
    QuitProgram();
  }

  w = (WindoW) MakeStreamFilteringWindow (input_file, output_file, is_binary, is_batch, is_submit);
  
  Show (w);
  Update ();

  ProcessEvents ();


  ArrowCursor ();
  Update ();


  return 0;
}


