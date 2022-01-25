/*   udvseq.c
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
* File Name:  udvseq.c
*
* Author:  Patrick Durand
*
* Version Creation Date:   5/3/99
*
* $Revision: 6.2 $
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

#include <udvseq.h>

/*******************************************************************************

  Function : WhatBspName()
  
  Purpose : retrieve the name of the bsp. Currently call only BioseqGetTitle;
  			may be changed in the future
  
  Parameters : bsp; pointer to the Bioseq
  				szName; filled with the Bioseq's name
  
  Return value : none

*******************************************************************************/
static void  WhatBspName(BioseqPtr bsp,CharPtr szName)
{
	szName=BioseqGetTitle(bsp);
}

/*******************************************************************************

  Function : WhatBspID()
  
  Purpose : retrieve the identifier of the bsp. 
  
  Parameters : bsp; pointer to the Bioseq
  				szAccess; filled with the Bioseq's accession number
  
  Return value : none

*******************************************************************************/
static void  WhatBspID(BioseqPtr bsp,CharPtr szAccess)
{
	SeqIdWrite(bsp->id,szAccess,PRINTID_TEXTID_ACCESSION,20);
}


/*******************************************************************************

  Function : WhatBspDataType()
  
  Purpose : retrieve the coding type of the sequence.
  
  Parameters : bsp; pointer to the Bioseq
  				szDtype; filled with the Bioseq's coding type
  
  Return value : none

*******************************************************************************/
static void  WhatBspDataType(BioseqPtr bsp,CharPtr szDtype)
{
Uint1 i;
Char *datatype[]={"IUPACna","IUPACaa","NCBI2na","NCBI4na",
					"NCBI8na","NCBIpna","NCBI8aa","NCBIeaa",
					"NCBIpaa","iupacaa3","NCBIstdaa","Unknown"};
	
	i=bsp->seq_data_type;

	switch(i){
		case Seq_code_iupacna:
		case Seq_code_iupacaa:
		case Seq_code_ncbi2na:
		case Seq_code_ncbi4na:
		case Seq_code_ncbi8na:
		case Seq_code_ncbipna:
		case Seq_code_ncbi8aa:
		case Seq_code_ncbieaa:
		case Seq_code_ncbipaa:
		case Seq_code_iupacaa3:
		case Seq_code_ncbistdaa:
			StringCpy(szDtype,datatype[i-1]);
			break;
		default:
			StringCpy(szDtype,datatype[11]);
			break;
	}
}

/*******************************************************************************

  Function : WhatBspStrand()
  
  Purpose : retrieve the strand def. of the BioSeq
  
  Parameters : bsp; pointer to the Bioseq
  				szStrand; filled with the Bioseq's strand def.
  
  Return value : none

*******************************************************************************/
static void  WhatBspStrand(BioseqPtr bsp,CharPtr szStrand)
{
Uint1 i;
Char *strand[]={"Single","Double","Mixed","Unknown"};
	
	i=bsp->strand;
	switch(i){
		case 1:
		case 2:
		case 3:
			StringCpy(szStrand,strand[i-1]);
			break;
		default:
			StringCpy(szStrand,strand[3]);
			break;
	}
}



/*******************************************************************************

  Function : WhatBspTopo()
  
  Purpose : retrieve the topology of the Bioseq
  
  Parameters : bsp; pointer to the Bioseq
  				szTopo; filled with the Bioseq's topology.
  
  Return value : none

*******************************************************************************/
static void  WhatBspTopo(BioseqPtr bsp,CharPtr szTopo)
{
Uint1 i;
Char *topo[]={"Linear","Circular","Tandem","Unknown"};
	
	i=bsp->topology;
	switch(i){
		case 1:
		case 2:
		case 3:
			StringCpy(szTopo,topo[i-1]);
			break;
		default:
			StringCpy(szTopo,topo[3]);
			break;
	}
}

/*******************************************************************************

  Function : WhatBspMol()
  
  Purpose : what's the bioseq ?
  
  Parameters : bsp; pointer to the Bioseq
  				szMol; filled with the Bioseq's type.
  
  Return value : TRUE if Bioseq is nucleic acid

*******************************************************************************/
static Boolean  WhatBspMol(BioseqPtr bsp,CharPtr szMol)
{
Uint1 i;
Char *mol[]={"DNA","RNA","AA","NA","Unknown"};				
	
	i=bsp->mol;
	switch(i){
		case Seq_mol_dna:
		case Seq_mol_rna:
		case Seq_mol_aa:
		case Seq_mol_na:
			StringCpy(szMol,mol[i-1]);
			break;
		default:
			StringCpy(szMol,mol[4]);
			break;
	}
	
	return(ISA_na(i));	
}

/*******************************************************************************

  Function : WhatBspRepr()
  
  Purpose : What's the bioseq representation
  
  Parameters : bsp; pointer to the Bioseq
  				szRepr; filled with the Bioseq's representation.
  
  Return value : none

*******************************************************************************/
static void  WhatBspRepr(BioseqPtr bsp,CharPtr szRepr)
{
Uint1 i;
Char *repr[]={"Virtual","Raw","Segmented","Constructed",
				"Reference","Consensus","Map","Delta","Unknown"};

	i=Bioseq_repr(bsp);
	switch(i){
		case Seq_repr_virtual:
		case Seq_repr_raw:
		case Seq_repr_seg:
		case Seq_repr_const:
		case Seq_repr_ref:
		case Seq_repr_consen:
		case Seq_repr_map:
		case Seq_repr_delta:
			StringCpy(szRepr,repr[i-1]);
			break;
		default:
			StringCpy(szRepr,repr[8]);
			break;
	} 
}

/*******************************************************************************

  Function : UDV_ReadBspDataForViewer()
  
  Purpose : call all of the above functions; exported function
  
  Parameters : 
  
  Return value : 

*******************************************************************************/
NLM_EXTERN void  UDV_ReadBspDataForViewer(BspInfoPtr bsp_i)
{
	WhatBspName(bsp_i->bsp,bsp_i->bspName);
	WhatBspID(bsp_i->bsp,bsp_i->bspAccNum);
	WhatBspDataType(bsp_i->bsp,bsp_i->bspDataType);
	WhatBspStrand(bsp_i->bsp,bsp_i->bspStrand);
	WhatBspTopo(bsp_i->bsp,bsp_i->bspTopo);
	bsp_i->bspMolNuc=WhatBspMol(bsp_i->bsp,bsp_i->bspMol);
	WhatBspRepr(bsp_i->bsp,bsp_i->bspRepr);
	bsp_i->bspLength=BioseqGetLen(bsp_i->bsp);
}

/*******************************************************************************

  Function : UDV_FreeListParaG()
  
  Purpose : delete the ParaG population
  
  Parameters : head of the ParaG val_node list
  
  Return value : none

*******************************************************************************/
NLM_EXTERN void UDV_FreeListParaG(ValNodePtr PNTR vnp_head)
{
ValNodePtr vnp;
ParaGPtr pgp;

	if ((*vnp_head)){
		for(vnp=(*vnp_head) ; vnp != NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if (pgp->pFeatList) ValNodeFree(pgp->pFeatList);
				Free(pgp);
			}		
		}	

		ValNodeFree(*vnp_head);
		*vnp_head=NULL;
	}
}

/*******************************************************************************

Function: UDV_DecodeIdxFeat()

Purpose: decode a 32 bits integer containing val1 (lower 16 bits)
		and val2 (higher 16 bits)

Parameters: index_g; value to decode in 'val1' and 'val2'

*******************************************************************************/

NLM_EXTERN void  UDV_DecodeIdxFeat (Uint4 index_g, Uint2Ptr val1,
									   Uint2Ptr val2)
{
Uint2Ptr  index_g2;

	index_g2 = (Uint2Ptr) (&index_g);
	*val1 = (Uint2) index_g2 [0];
	*val2 = (Uint2) index_g2 [1];
}

/*******************************************************************************

Function: UDV_EncodeIdxFeat()

Purpose: encode a 32 bits integer with val1 (lower 16 bits)
		and val2 (higher 16 bits)

Parameters: values to decode -> 'val1' and 'val2'

Return : the 32 bits encoded value

*******************************************************************************/

NLM_EXTERN Uint4  UDV_EncodeIdxFeat (Uint2 val1,Uint2 val2)
{
Uint2 index_g[2];
	
	index_g[0]=val1;
	index_g[1]=val2;
	
	return *((Int4Ptr) index_g);
	
}

/*******************************************************************************

Function: UDV_IsTranslationNeeded()

Purpose: for a CDS in a ParaG, determine whether a translation is needed or not

Parameters: context ; feature data
			pgp ; ParaG structure

Return : TRUE is a translation is needed

*******************************************************************************/
NLM_EXTERN Boolean UDV_IsTranslationNeeded(SeqMgrFeatContextPtr context,ParaGPtr pgp)
{
Int2 i,numivals2,i_decal;	/*counters*/
Uint1 strand;

	strand=context->strand;
	
	/*temporary situation; will be modified in the future*/
	if (strand>Seq_strand_minus ||
		strand==Seq_strand_unknown) strand=Seq_strand_plus;

	/*strand PLUS*/
	if (strand==Seq_strand_plus){
		numivals2=context->numivals*2;
		i=0;
		i_decal=2;
	}

	/*strand MINUS*/
	if (strand==Seq_strand_minus){
		numivals2=2*context->numivals-2;
		i=numivals2;
		i_decal=-2;
	}		

	while (TRUE){
		/*if ivals.stop<= start ParaG : not yet in the current ParaG*/
		if (context->ivals[i+1]<pgp->StartLetter) {
			if (strand==Seq_strand_plus){
				if (numivals2>2 && i+2<numivals2){
					/*stop ParaG < start next ivals -> inter-region: fill 
					the ParaG with a thin line; this is the case
					for coding region: draw thin line to delineate the introns*/		
					if (context->ivals[i+2]>pgp->StopLetter){
						return(FALSE);
					}
				}
			}
			if (strand==Seq_strand_minus){
				if (numivals2>2 && i-2>-1){
					/*stop ParaG < start next ivals -> inter-region: fill 
					the ParaG with a thin line; this is the case
					for coding region: draw thin line to delineate the introns*/		
					if (context->ivals[i-2]>pgp->StopLetter){
						return(FALSE);
					}
				}
			}
		}

		if (strand==Seq_strand_plus){
			i=i+i_decal;
			if (i>numivals2-2) break;
		}
		if (strand==Seq_strand_minus){
			i=i+i_decal;
			if (i<0) break;
		}
	}		
	return(TRUE);
}



/*******************************************************************************

Function: UDV_ParaGFTableFeatures()

Purpose: retrieve features for several ParaG at a time

Parameters: see explore.h

*******************************************************************************/

NLM_EXTERN Boolean LIBCALLBACK UDV_ParaGFTableFeatures (SeqFeatPtr sfp, 
			SeqMgrFeatContextPtr context)

{
ParaGFeaturesInLocPtr pgfl;
ValNodePtr vnp,vnp2,vnp3,vnp4;
ParaGPtr pgp;
Int4 BeginLines=0;
Int2 nLines=0,i=1;
Boolean FeatInParaG=FALSE;
Boolean bFirst=TRUE,bTrouve=FALSE;
Boolean IsTransNeeded=FALSE;

	if (!context->sfp) return (TRUE);

	pgfl = (ParaGFeaturesInLocPtr) context->userdata;	

	vnp=pgfl->ParaG_head;
	pgp=(ParaGPtr)vnp->data.ptrvalue;

	BeginLines=pgp->StartLine;

	pgfl->nFeat++;
	
	/*HET feature correction for the ends*/
	if(context->featdeftype== FEATDEF_HET){
		context->right=context->ivals[2*context->numivals-1];
	}
	
	while(TRUE){		
		pgp=(ParaGPtr)vnp->data.ptrvalue;
		
		FeatInParaG=FALSE;
		if (pgp){
			if (pgp->StartLetter>context->right){
				pgfl->ParaG_last_head=vnp;/*uses to count last ParaG of the list
											where there is no feature*/
				break;
			}
			/*Feature in pgp ?*/
			if(pgp->StopLetter>=context->left && pgp->StopLetter<=context->right){
					FeatInParaG=TRUE;
			}
			else if(pgp->StartLetter>=context->left && 
				pgp->StartLetter<=context->right){
					FeatInParaG=TRUE;
			}
			else if(pgp->StartLetter<=context->left && 
				pgp->StopLetter>=context->right){
					FeatInParaG=TRUE;

			}

			if (FeatInParaG){
				if (bFirst) {
					pgfl->ParaG_next_head=vnp;
					bFirst=FALSE;
				}
				nLines=1;
				if (context->sfp->data.choice==SEQFEAT_CDREGION) {
					IsTransNeeded=UDV_IsTranslationNeeded(context,pgp);
					if (IsTransNeeded) {
						nLines++;
						pgp->nTrans++;
					}
				}
				if (pgfl->ShowFeatures) {
					if (context->left<=pgp->OccupyTo){/*new line of features*/
						pgp->nLines+=nLines;
						pgp->nFeatLines+=nLines;
					}
					else{/*same feature line, but CDS*/
						/*for the translation*/
						if (context->sfp->data.choice==SEQFEAT_CDREGION){
							pgp->nLines++;
							pgp->nFeatLines++;
						}
					}
					pgp->OccupyTo=context->right;
				}
				pgp->nFeat++;
				/*(re)populate Features structure; if a previous structure already
				exists, use it; avoid a lot of MemNew/MemFree; see also the function
				which create/(re)populate ParaG structure (CreateParaGList)*/
				if (pgp->pFeatList==NULL) {
					pgp->pFeatList = ValNodeAddInt(&pgp->pFeatList,1,
							UDV_EncodeIdxFeat ((Uint2) context->itemID,
							(Uint2) context->index));
					if (!pgp->pFeatList) return TRUE;
				}
				else {
					i=1;bTrouve=FALSE;
					for (vnp2=pgp->pFeatList ; vnp2!=NULL ; vnp2=vnp2->next){
						if (i==pgp->nFeat){
							bTrouve=TRUE;
							break;
						}
						i++;
						vnp3=vnp2;
					}
					/*mem alloc only if it's necessary*/
					if (bTrouve) {
						vnp2->data.intvalue=UDV_EncodeIdxFeat((Uint2)
								context->itemID,(Uint2) context->index);
					}
					else{
						vnp4=ValNodeAddInt(&vnp3,1,UDV_EncodeIdxFeat((Uint2)
								context->itemID,(Uint2) context->index));
						if (!vnp4)return TRUE;
					}
				}
			}

			pgp->StartLine=BeginLines;

			BeginLines+=pgp->nLines;

			vnp=vnp->next;

			if (vnp==NULL){
				pgfl->ParaG_next_head=NULL; /*end*/
				break;
			}
		}
	}

	pgfl->nTotLines_new=BeginLines;
	
	if (pgfl->ParaG_next_head)
		pgfl->ParaG_head=pgfl->ParaG_next_head;

	return TRUE;
}


/*******************************************************************************

Function: UDV_CreateOneFeatureIndex()

Purpose: Call Sequence Manager to create the Feature Index of a BioSeq.

Parameter: entityID_seq ; entityID of the BioSeq.

Return value: the index value to use with others Sequence Manager functions.
		(Ex.: SeqMgrExplore... functions)

*******************************************************************************/

NLM_EXTERN Uint2  UDV_CreateOneFeatureIndex(Uint2 entityID_seq, BioseqPtr bsp)
{
Uint2 EntityID;

	EntityID=entityID_seq;
	/*Is a previous Index exists for this sequence*/
	if (SeqMgrFeaturesAreIndexed (entityID_seq) == 0) {
		/*create the index*/
		if (bsp)
			EntityID=SeqMgrIndexFeatures (0, bsp);
		else
			EntityID=SeqMgrIndexFeatures (entityID_seq, NULL);
	}
	
	if (EntityID == 0) {
		Message (MSG_ERROR, "SeqMgrIndexFeatures failed.");
		EntityID=INDEX_CREATION_ERROR;
	}
	return(EntityID);
}

/*******************************************************************************

  Function : UDV_CreateParaGList()
  
  Purpose : Create and populate ParaG val_node_list without features
  
  Parameters : 	nCharByLine; # of char by line
  				bsp_length; length of the bsp
				ShowTop; show scale on top if TRUE
				ShowTick; show scale's ticks if TRUE
				ShowSequence; show sequence if TRUE
				ShowBlank; add blank line at the bottom of the ParaG
				nTotL; total number of single lines for the whole ParaG list
				ParaG_head; head of the ParaG list

  Note : use this function to create a simple representation of a sequence
  				without the features. If you want to show the features, first
				call this function, then UDV_CreateOneFeatureIndex() to create
				feature index, and finally 	PopulateParaGFeatures() to populate
				ParaG with the features.		

  Return value : a ValNodePtr which is the head of the ParaG list
			Indeed, in this function ,the value of ParaG_head (passed in here) is
			modified. This is done to speed up the population of ParaG.
			
			(see also nTotL)
			
*******************************************************************************/
NLM_EXTERN ValNodePtr UDV_CreateParaGList(Int2 nCharByLine,
			Int4 bsp_length,Int4 from,Int4 to,
			Boolean ShowTop,Boolean ShowTick,Boolean ShowSequence, 
			Boolean ShowBlank,Int4Ptr nTotL,ValNodePtr ParaG_head)
{
ParaGPtr 	pgp=NULL;				/*ParaG data*/
ValNodePtr 	local_head=NULL,		/*these are used to scan ParaG*/
			vnp1=NULL,
			vnp2=NULL;				
Boolean		isFailed=FALSE;			/*memory allocation failure*/
Boolean		bFirst=TRUE;			/*used to avoid non-neces. mem alloc*/
Int4		start=0,				/*letter start, zero based*/
			stop=0,					/*letter stop. zero based*/
			nTotLines=0;			/*Total lines in viewer*/
Int1		minLineByParaG=0;		/*min height of a paraG, nb. of lines*/			
Int4		nParaG=0;				/*ParaG counter*/
			
	/*Minimum ParaG height*/
	if (ShowTop) minLineByParaG++;
	if (ShowTick) minLineByParaG++;
	if (ShowSequence) minLineByParaG++;
	if (ShowBlank) minLineByParaG++;

	/*stop=MIN(nCharByLine,bsp_length);*/

	local_head=ParaG_head;
	vnp1=ParaG_head;
	start=from;
	/*to++;
	stop=MIN(nCharByLine,to);*/
	
	/*(re)populate ParaG structure; if a previous structure already
	exists, use it; avoid a lot of MemNew/MemFree*/
	while(!(start>=to)){
		if (vnp1){
			pgp=(ParaGPtr)vnp1->data.ptrvalue;

			if (!pgp){
				isFailed=TRUE;
				break;
			}

			/*Features; note that pgp->FeatList is not deleted : see the
			function which populates this structure (ParaGFTableFeatures)*/
			pgp->nFeat=0;pgp->nTrans=0;
			if (bFirst) bFirst=FALSE;
			vnp2=vnp1;
			vnp1=vnp1->next;
		}
		else{
			/*allocate new ParaG*/
			pgp=(ParaGPtr)MemNew(sizeof(ParaG));
			if (!pgp){
				isFailed=TRUE;
				break;
			}
			MemSet(pgp,0,sizeof(ParaG));
			
			/*create a new node*/
			if (vnp2==NULL) vnp2=ValNodeAddPointer(NULL,0,pgp);
			else vnp2=ValNodeAddPointer(&vnp2,0,pgp);

			if (!vnp2){
				isFailed=TRUE;
				break;
			}
			
			if (bFirst){
				local_head=vnp2;
				bFirst=FALSE;
			}
		}

		/*Fill in pgp*/
		pgp->StartLine=nTotLines;
		pgp->nLines=minLineByParaG;
		pgp->StartLetter=start;
		pgp->NumOrder=nParaG+1;
		pgp->nFeatLines=0;
		
		/*modify values*/
		nTotLines+=minLineByParaG;
		start+=nCharByLine;
		stop=start-1;/*nCharByLine;*/
		
		if (stop>to/*bsp_length*/) stop=to/*bsp_length*/;
		pgp->StopLetter=stop;
		pgp->OccupyTo=pgp->StopLetter;

		nParaG++;
	}
	
	if (isFailed){
		UDV_FreeListParaG(&local_head);
		nTotLines=0;		
	}

	if (vnp1){
		ValNodeFree(vnp1);
		vnp2->next=NULL;
	}

	*nTotL=nTotLines;

	return(local_head);
}

/*******************************************************************************

  Function : UDV_PopulateParaGFeatures()
  
  Purpose : populate ParaG val_node_list with features
  
  Parameters : 	GrData; graphical data
  				bsp; BioseqPtr of the Bioseq
				ParaG_vnp; valnode list of the ParaG. Each node has a
					data.ptrvalue which is of type ParaGPtr.
				rcP; size and position of the UnD viewer's panel
				nTotL; total number of single lines for the whole ParaG list

  Note : this function, if called, MUST BE used after a call to	
  		UDV_CreateParaGList().	

  Return value : TRUE if features found; otherwise FALSE (see also nTotL)

*******************************************************************************/
NLM_EXTERN Boolean UDV_PopulateParaGFeatures(BioseqPtr bsp,ValNodePtr ParaG_vnp,
			Boolean ShowFeatures,Int4Ptr nTotL)
{
ParaGFeaturesInLoc pgfl;	/*used for the Explore features function*/
ValNodePtr vnp;				/*to scan ParaG list*/
ParaGPtr pgp;				/*to modify ParaG, when needed*/
Int4 n=0;					/*a little counter*/
Boolean nRet;

	if (ParaG_vnp == NULL) {*nTotL=0;return(FALSE);}

	/*prepare the Exploration of Features*/
	pgfl.ParaG_head=ParaG_vnp;
	pgfl.nTotLines_new=0;
	pgfl.ParaG_next_head=NULL;
	/*pgfl.LineH=GrData.udv_font.LineHeight;*/
	/*pgfl.rcP_top=rcP.top;*/
	pgfl.ParaG_last_head=NULL;
	pgfl.ShowFeatures=ShowFeatures;
	
	/*Explore the features */
	SeqMgrExploreFeatures (bsp, (Pointer) &pgfl, UDV_ParaGFTableFeatures, 
				NULL, NULL, NULL);

	/*complete the list where no features have been found; actually, this
	code is used to locate (RecT) correctly the end of the ParaG list*/
	if (pgfl.ParaG_next_head && pgfl.ParaG_last_head){
		n=pgfl.nTotLines_new;
		for(vnp=pgfl.ParaG_last_head ; vnp != NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				pgp->StartLine=n;
				n+=pgp->nLines;
			}		
		}
		pgfl.nTotLines_new=n;	
	}
	
	/*pgfl.nTotLines_new could be 0 if no features*/
	if (pgfl.nTotLines_new>0) {
		*nTotL=pgfl.nTotLines_new;
		nRet=TRUE;
	}
	else nRet=FALSE;

	return(nRet);
}

/*******************************************************************************

  Function : UDV_Read_Sequence()
  
  Purpose : read a sequence 
  
  Parameters : Warning, from must be < than to; always
  
  Return value : the sequence 

*******************************************************************************/
NLM_EXTERN CharPtr UDV_Read_Sequence (SeqIdPtr sip, Int4 from, Int4 to, 
		Boolean IsProt,Int2 len)
{
/*BioseqPtr        bsp;*/
SeqLocPtr        slp;
SeqPortPtr       spp;
CharPtr          str = NULL;
Uint1			residue;
Uint2			i=0;

	/*from always < than to*/
	slp = SeqLocIntNew (from, to, Seq_strand_plus, sip);
	spp = SeqPortNewByLoc (slp, (Uint1)(IsProt==TRUE ? Seq_code_iupacaa 
			: Seq_code_iupacna));
	if (spp != NULL) {
		str = MemNew ((len+1) * sizeof(Char));
		if (!str) return(NULL);
		while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF) {
			if (IS_residue(residue)) {
				str[i] = residue;
				i++;
			}
		}
		SeqPortRead(spp, (BytePtr)str, len);
		SeqPortFree (spp);
	}   

	SeqLocFree (slp);
	
	return str;
}
