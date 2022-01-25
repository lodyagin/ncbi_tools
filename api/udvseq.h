/*   udvseq.h
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
* File Name:  udvseq.h
*
* Author:  Patrick Durand
*
* Version Creation Date:   5/3/99
*
* $Revision: 6.1 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#ifndef _UDVSEQ_
#define _UDVSEQ_

#ifdef __cplusplus
extern "C" {
#endif

#include <explore.h>
#include <ncbi.h>
#include <objfdef.h>
#include <objloc.h>
#include <objseq.h>
#include <seqport.h>
#include <sequtil.h>

/*******************************************************************************

  The major element of this module is: paragraph (or ParaG)

  This module is Vibrant free and can be used by non-vibrant software.
  
*******************************************************************************/

	/*****************************************************************************

	Structure of a ParaG
	
              10         20        30   Line 0
	      .    |     .    |    .    |   Line 1
	1 GFDEDJSHHF DHGCJFDHCGNJFGJDHFGJ   Line 3
        |===>......===>.....===>|       Line 4
		/\/\/\/\    >------<  ^         Line n
		
		Line 0: numerical scale
		Line 1: ticks (minor: . ; major: |)
		Line 3: left scale (if needed) + the sequence
		Line 4: features
		Line n: features
	*****************************************************************************/


/*******************************************************************************

	DEFINES

*******************************************************************************/

/*Error value; if Feature Index failed*/
#define INDEX_CREATION_ERROR FALSE

/*******************************************************************************

	STRUCTURES

*******************************************************************************/

typedef struct parag {/*Paragraph information*/
	Int4 NumOrder;
	/*ParaG graphical values*/
	Int4 StartLine;				/*this ParaG starts at this line*/
	Int4 nLines;				/*and contains nLines (scale+seq+feat)*/
	Int2 nFeatLines;			/*number of lines with features*/
	Int4 StartLetter;			/*first letter of the bioseq to show*/
	Int4 StopLetter;			/*last letter of the bioseq to show*/

	/*ParaG Features List*/
	ValNodePtr pFeatList;		/*list of itemID,index (Feature Index
								values)*/
	Int2 nFeat;					/*number of features*/
	Int2 nTrans;				/*number of translation*/
	Int4 OccupyTo;				/*used to populate*/
	} ParaG, PNTR ParaGPtr;

/*structure used to initialize ParaG with features*/
typedef struct paragfeaturesinloc{
	ValNodePtr  ParaG_head;
	Int4		nTotLines_new;
	ValNodePtr	ParaG_next_head;
	ValNodePtr	ParaG_last_head;
	Int2		LineH;
	/*Int2		rcP_top;*/
	Int4		nFeat;
	Boolean 	ShowFeatures;
	} ParaGFeaturesInLoc,PNTR ParaGFeaturesInLocPtr;

typedef struct bspinfo {
		/*general data*/
	Uint2			bsp_entityID;
	Uint2			bsp_itemID;
	Uint2			bsp_itemType;
	BioseqPtr 		bsp;
	Char			bspName[41];
	Char			bspAccNum[21];
	Char			bspRepr[21];
	Char			bspMol[21];
	Boolean         bspMolNuc;
	Char			bspTopo[21];
	Char 			bspStrand[21];
	Char 			bspDataType[21];
	Int4			bspLength;
		/*Sequence Buffer*/
	Int4 			StartBuf;	/*buffer start here in the sequence*/
	Int4 			StopBuf;	/*buffer stop here... (0 based values)*/
	Int2 			LengthBuf;	/*buffer size*/
	CharPtr 		SeqBuf;		/*buffer sequence*/
	ValNodePtr 		PgpStartBuf;/*Pgp where start buffer is located*/
	}BspInfo, PNTR BspInfoPtr;

/*******************************************************************************

	FUNCTIONS DECLARATION
	(see .c module for a complete description)

*******************************************************************************/

/*BSP information*/
NLM_EXTERN void  UDV_ReadBspDataForViewer(BspInfoPtr bsp_i);
/*Feature management*/
NLM_EXTERN void UDV_DecodeIdxFeat (Uint4 index_g, Uint2Ptr val1,
		Uint2Ptr val2);
NLM_EXTERN Uint4  UDV_EncodeIdxFeat (Uint2 val1,Uint2 val2);
NLM_EXTERN Boolean UDV_IsTranslationNeeded(SeqMgrFeatContextPtr context,
		ParaGPtr pgp);
NLM_EXTERN Boolean LIBCALLBACK UDV_ParaGFTableFeatures (SeqFeatPtr sfp, 
		SeqMgrFeatContextPtr context);
NLM_EXTERN Uint2 UDV_CreateOneFeatureIndex(Uint2 entityID_seq, 
		BioseqPtr bsp);
NLM_EXTERN void UDV_FreeListParaG(ValNodePtr PNTR vnp_head);
NLM_EXTERN ValNodePtr UDV_CreateParaGList(Int2 nCharByLine,
		Int4 bsp_length,Int4 from,Int4 to,
		Boolean ShowTop,Boolean ShowTick,Boolean ShowSequence, 
		Boolean ShowBlank,Int4Ptr nTotL,ValNodePtr ParaG_head);
NLM_EXTERN Boolean UDV_PopulateParaGFeatures(BioseqPtr bsp,
		ValNodePtr ParaG_vnp,Boolean ShowFeatures,Int4Ptr nTotL);
/*Sequence reader*/
NLM_EXTERN CharPtr UDV_Read_Sequence (SeqIdPtr sip, Int4 from, Int4 to, 
		Boolean IsProt,Int2 len);



#ifdef __cplusplus
}
#endif

#endif /* ndef _UDVSEQ_ */

