/* vastchkid.c
 *
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
 * File Name: vastchkid.c
 *
 * Author: Yanli Wang 
 *
 * Version Creation Date: 6/16/98
 *
 * $Log: vastchkid.c,v $
 * Revision 6.1  1998/07/17 18:42:33  madej
 * Handles local sequence IDs.
 *
 *
 */
#include <stdio.h>
#include <ncbi.h>
#include <accentr.h>
#include <objalign.h>
#include <objseq.h>
#include <objmgr.h>
#include <lsqfetch.h>
#include <netentr.h>
#include <www.h>
#include <sys/resource.h>
#include <mmdbapi.h>
#include <mmdbapi1.h>
#include "vastlocl.h"
#include "mmdblocl.h"
#include "mmdbdata.h"
#include "vast2mage.h"
#include "vastsrv.h"
#include <asnmime.h>
#include <objmime.h>
#include <strimprt.h>

#define Local_max 20 
Char pSegmentMaster[Local_max];
/*--------------- find PDBSeqId ------------------*/
PDBSeqIdPtr GetPdbSeqId(SeqIdPtr sip)
              /* may need to be modified according to how bioseq id is */
{
  SeqIdPtr  seq_id = NULL;
  PDBSeqIdPtr pdb_seq_id = NULL;

  seq_id = sip;

  while(seq_id != NULL){
     if(seq_id->choice == 15){
        pdb_seq_id = seq_id->data.ptrvalue;
        break;
     }
     seq_id = seq_id->next;
  }

  return(pdb_seq_id);
}
/*---------------- remove redundancy sequences ---------------*/ 
SeqEntryPtr ExtractThisSep(SeqEntryPtr sepThis, CharPtr PDBName, Char Chain)
{
  SeqEntryPtr sepThisNew = NULL;
  BioseqPtr bsp = NULL;
  BioseqSetPtr bssp = NULL;
  PDBSeqIdPtr pdb_seq_id;

  Int2 choice;
  Char ThisChain;

  sepThisNew = SeqEntryNew();
  sepThisNew->choice = 1;

  choice = sepThis->choice;
  if(choice == 1) {
     sepThisNew->data.ptrvalue = sepThis->data.ptrvalue;
     return (sepThisNew);
  }
  else if(choice == 2){
     bssp = sepThis->data.ptrvalue;
     sepThis = bssp->seq_set;

     while(sepThis != NULL){
        bsp = sepThis->data.ptrvalue;
        pdb_seq_id = GetPdbSeqId(bsp->id);
        if(pdb_seq_id->chain !=NULL) ThisChain = (Char) pdb_seq_id->chain;
        else ThisChain = ' ';

        if(ThisChain == Chain) {
           sepThisNew->data.ptrvalue = bsp;
           return (sepThisNew);
        }

        sepThis = sepThis->next;
     }
  }

}
/*---------------- make a new SeqEntry ------------------------*/
NcbiMimeAsn1Ptr ScreenSequence( NcbiMimeAsn1Ptr pvnNcbi, CharPtr JobID)
{
  BiostrucAlignPtr pbsaStruct = NULL;
  BiostrucAnnotSetPtr pbsaThis = NULL, pbsaHead = NULL;
  BiostrucFeatureSetPtr pbsfsThis = NULL;
  BiostrucFeaturePtr pbsfThis = NULL;
  SeqEntryPtr sepHead = NULL, sepHeadNew = NULL, sepThis = NULL;

  Char PDBName[Local_max], Chain, Domain;
  Char pSegmentSlave[Local_max];

  pbsaStruct = pvnNcbi->data.ptrvalue;

  pbsaThis = pbsaStruct->alignments;
  pbsfsThis = pbsaThis->features;
  pbsfThis = pbsfsThis->features;

  StringNCpy(pSegmentMaster, pbsfThis->name, 6);
  StringNCpy(PDBName, pSegmentMaster, 4);
  Chain = pSegmentMaster[4];
  
  sepHead = pbsaStruct->sequences;
  if(JobID){
     sepThis = SeqEntryNew();
     sepThis->choice = 1;
     sepThis->data.ptrvalue = sepHead->data.ptrvalue;
   }
  else sepThis = ExtractThisSep(sepHead, PDBName, Chain); 
  ValNodeLink(&sepHeadNew, sepThis);

  sepHead = sepHead->next;
  while(pbsfThis){
     StringCpy(pSegmentSlave, pbsfThis->name + 7);
     pSegmentSlave[6]='\0';

     StringNCpy(PDBName, pSegmentSlave, 4);
     Chain = pSegmentSlave[4];

     sepThis = ExtractThisSep(sepHead, PDBName, Chain);
     ValNodeLink(&sepHeadNew, sepThis);

     pbsfThis = pbsfThis->next;
     sepHead = sepHead->next; 
  }

  pbsaStruct->sequences = sepHeadNew; 

  return(pvnNcbi);
}
/*---------------- replace original seq_id by local_id  --------------*/
void AssignLocaIdToBiostruc(SeqIdPtr sip, Char *PDBName, Char Chain,BiostrucPtr bsp)
{
  BiostrucDescrPtr pbsdrThis = NULL;
  MoleculeGraphPtr pmgThis = NULL;
  BiomolDescrPtr pbmdrThis = NULL;

  Char ThisPDBName[Local_max], *ThisChain;

  StringCpy(ThisPDBName, PDBNAME_DEFAULT);
  pbsdrThis = ValNodeFindNext(bsp->descr, NULL, BiostrucDescr_name);
  if(pbsdrThis) StringNCpy(ThisPDBName, pbsdrThis ->data.ptrvalue, 4);

/*if(StringCmp(PDBName, ThisPDBName) != 0) {
     printf("Something is wrong!\n");
     return;
  }  */

  pmgThis = bsp->chemical_graph->molecule_graphs;

  while(pmgThis) {
                     /* loop over molecule_graphs to match chain */ 
     pbmdrThis = ValNodeFindNext(pmgThis->descr, NULL, BiomolDescr_name);
     if(pbmdrThis){
        ThisChain = StringSave((CharPtr)pbmdrThis->data.ptrvalue); 
     }
     else ThisChain = StringSave(STRING_UNK_MOL);
     if(Chain == *ThisChain) {
        pmgThis->seq_id = sip; 
        return;
      }

     pmgThis = pmgThis->next;
  }

}
/***************************************************************************/                                                                       
/*            GetObjectId() Ken                                            */
/***************************************************************************/
ObjectIdPtr GetObjectId(BiostrucPtr bsp, CharPtr pSegmentThis)
{
  ObjectIdPtr oip;
  SeqIdPtr sip;
  MoleculeGraphPtr mgp;
  
  for (mgp = bsp->chemical_graph->molecule_graphs; mgp; mgp = mgp->next)
  {
    sip = mgp->seq_id;
    oip = (ObjectIdPtr)sip->data.ptrvalue;
    if (oip->str[5] == pSegmentThis[4]) break;
  }
  
  return oip;
}
/*--------------- check SeqId for each Bioseq -----------------*/
Boolean ReplaceBioseqId(SeqIdPtr sip, Char *PDBName, Char Chain, BioseqPtr bsp)
{
  PDBSeqIdPtr pdb_seq_id;
  Char ThisChain;
  
  pdb_seq_id = GetPdbSeqId(bsp->id);
          /* may need to be modified according to how bioseq id is */
          /* in bioseq fetched from Entrez, PDBSeqId exists, but not sure */
          /* how struc2seq works */
  if(pdb_seq_id->chain !=NULL) ThisChain = (Char) pdb_seq_id->chain;
  else ThisChain = ' ';
  
  if(ThisChain == Chain) {
     bsp->id = sip;
     return TRUE;
  }

  else return FALSE;

}
/*---------------- replace original seq_id by local_id  --------------*/
void AssignLocaIdToBioseq(SeqIdPtr sip, Char *PDBName, Char Chain, SeqEntryPtr sep)
{
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  Int2 choice;
 
  Boolean IdReplaced = FALSE;

  choice = sep->choice;
  if(choice == 1) {
     bsp = sep->data.ptrvalue;
     ReplaceBioseqId(sip, PDBName, Chain, bsp);
   }
   else if(choice == 2){
      bssp = sep->data.ptrvalue;
      sep = bssp->seq_set;

      while(sep != NULL){
         bsp = sep->data.ptrvalue;
         IdReplaced = ReplaceBioseqId(sip, PDBName, Chain, bsp); 
         if(IdReplaced) return;
         sep = sep->next;
      }
   }
   
}
/*---------------- replace original seq_id by local_id  --------------*/
void AssignLocaIdToSeqAlign(SeqIdPtr sip, Char *PDBName, Char Chain, SeqAlignPtr salpHead)
{

  DenseDiagPtr ddp = NULL;
  SeqIdPtr  sipThis = NULL;
  PDBSeqIdPtr pdb_seq_id;
  Char ThisChain;
 
  ddp = salpHead->segs;

  while(ddp){
     sipThis = ddp->id;
     sipThis = sipThis->next;
              /* get to the slave SeqId */

     pdb_seq_id = GetPdbSeqId(sipThis);
     if(pdb_seq_id != NULL){
        if(pdb_seq_id->chain !=NULL) ThisChain = (Char) pdb_seq_id->chain;
        else ThisChain = ' ';
        if(ThisChain == Chain) ddp->id->next = sip;
                       /* actually match chain here may not be neccessary */
     }

     ddp = ddp->next; 

  }

}
/*---------------- To put Local_id in segs for Vastsearch - Ken -----------*/
void AssignLocaIdToSeqAlignMaster(SeqIdPtr sip, SeqAlignPtr salpHead)
{

  DenseDiagPtr ddp = NULL;
  SeqIdPtr  sipTemp = NULL;
  PDBSeqIdPtr pdb_seq_id;
  Char ThisChain;
 
  ddp = salpHead->segs;

  while(ddp)
  {
     sipTemp = ddp->id->next;
     ddp->id = sip;
     sip->next = sipTemp;
     ddp = ddp->next; 
  }

}
/*---------------- to check if an identical chain exist ----------------*/
Boolean FindIdSlave(Char *pSegmentThis, BiostrucFeaturePtr pbsfHead, Int4 iCount1)
{
  Int4 iCount2 = 0;

  if(StringNCmp(pSegmentThis, pSegmentMaster, 5) == 0){
                         /* only compare PDB code and Chain Id */
     return TRUE;
  }
  else if(iCount1 == 0) {
     return FALSE;
  }

  iCount2 = 0;
  while(pbsfHead){
     if(iCount2 >= iCount1) return FALSE;
     if(StringNCmp(pSegmentThis, pbsfHead->name + 7, 5) == 0) {
                          /* only compare PDB code and Chain Id */
        return TRUE;
     }

     iCount2++;
     pbsfHead = pbsfHead->next;
  }

  return FALSE;

} 
/*---------------- make new LocalId ------------------------*/
SeqIdPtr MakeLocalId(Char *str)
{
  SeqIdPtr sip;
  ObjectIdPtr oip;

  oip = ObjectIdNew();
  oip->str = str;

  sip = ValNodeNew(NULL);
  sip->choice = SEQID_LOCAL;
  sip->data.ptrvalue = oip;     

  return(sip);
}
/*---------------- to assign local id to identical chain ---------------*/
NcbiMimeAsn1Ptr CheckId(NcbiMimeAsn1Ptr pvnNcbi, CharPtr JobID)
{
            /* to assign local id to identical chain so that to work around */
            /* object manager for messaging */  /* yanli */

  BiostrucAlignPtr pbsaStruct = NULL;
  BiostrucPtr pbsMaster = NULL, pbsSlaveHead = NULL;
  BiostrucAnnotSetPtr pbsaThis = NULL, pbsaHead = NULL;
  SeqAnnotPtr psaAlignHead = NULL; 
  SeqEntryPtr sepHead;

  BiostrucFeatureSetPtr pbsfsThis = NULL, pbsfsHead = NULL;
  BiostrucFeaturePtr pbsfThis = NULL, pbsfHead = NULL;
  SeqAlignPtr salpHead;

  SeqIdPtr sip, sip_cpy, sipMaster, sip_temp;
  ObjectIdPtr oipMaster, oip;

  Int4 iCount1 = 0, iCount_ID = 0;
  Boolean IdChainFound = FALSE;

  Char pSegmentSlave[Local_max];
  Char PDBName[Local_max], Chain, Domain; 

  Char AlignIdStr[Local_max];
  Int4 AlignId;    

  typedef struct local_string{
  Char str[Local_max];
  }Local_String, PNTR Local_StringPtr;

  Local_String **str; 
 
  pvnNcbi = ScreenSequence(pvnNcbi, JobID);

  pbsaStruct = pvnNcbi->data.ptrvalue;
  
  pbsaThis = pbsaStruct->alignments; 
  pbsSlaveHead = pbsaStruct->slaves;
  sepHead = pbsaStruct->sequences;
  psaAlignHead = pbsaStruct->seqalign;

  pbsfsThis = pbsaThis->features;
  pbsfThis = pbsfsThis->features;
                     /* from vastsrv page, we are always considering one */
                     /* domain of the master, so no loop for pbsfsThis */
 
  StringNCpy(pSegmentMaster, pbsfThis->name, 6);
  pSegmentMaster[6]='\0';
/*pSegmentMaster[0]=pbsfThis->name[0];
  pSegmentMaster[1]=pbsfThis->name[1];
  pSegmentMaster[2]=pbsfThis->name[2];
  pSegmentMaster[3]=pbsfThis->name[3];
  pSegmentMaster[4]=pbsfThis->name[4];
  pSegmentMaster[5]=pbsfThis->name[5];  
  pSegmentMaster[6]='\0';  */
 
  /* see if Master Sequence was generated by VastSearch - Ken*/ 
  if (JobID)
  {
    for (salpHead = psaAlignHead->data; salpHead; salpHead = salpHead->next)
    {
      oipMaster = ObjectIdNew();
      oipMaster = GetObjectId(pbsaStruct->master, pSegmentMaster);
    
      sipMaster = ValNodeNew(NULL);
      sipMaster->choice = SEQID_LOCAL;
      sipMaster->data.ptrvalue = oipMaster;
      
      sip_temp = AsnIoMemCopy(sipMaster, (AsnReadFunc) SeqIdAsnRead, (AsnWriteFunc) SeqIdAsnWrite);
      AssignLocaIdToSeqAlignMaster(sip_temp, salpHead);
    }
  }
    
  salpHead = psaAlignHead->data;
             /* from SeqAnnot to SeqAlign */ 
  sepHead = sepHead->next;
              /* go to the slave sequences */ 

  iCount1 = 0;
  while(pbsfThis){
     StringCpy(pSegmentSlave, pbsfThis->name + 7);
     pSegmentSlave[6]='\0';

     IdChainFound = FindIdSlave(pSegmentSlave, pbsaThis->features->features, iCount1);
     if(IdChainFound) iCount_ID++;
     iCount1++;
     pbsfThis = pbsfThis->next;
  }

  str = (Pointer) MemNew((size_t) ((iCount_ID + 1) * sizeof(Pointer)));
  for(iCount1 = 0; iCount1 < iCount_ID; iCount1++){
     str[iCount1] = (Local_StringPtr) MemNew(sizeof(Local_String));
  }

  pbsfThis = pbsfsThis->features;
  iCount_ID = 0;
  iCount1 = 0;
  while(pbsfThis){

     StringCpy(pSegmentSlave, pbsfThis->name + 7);
     pSegmentSlave[6]='\0';

/*   pSegmentSlave[0]=pbsfThis->name[7];
     pSegmentSlave[1]=pbsfThis->name[8];
     pSegmentSlave[2]=pbsfThis->name[9];
     pSegmentSlave[3]=pbsfThis->name[10];
     pSegmentSlave[4]=pbsfThis->name[11];
     pSegmentSlave[5]=pbsfThis->name[12];     
     pSegmentSlave[6]='\0';   */

     IdChainFound = FindIdSlave(pSegmentSlave, pbsaThis->features->features, iCount1);
     if(IdChainFound){
        StringNCpy(str[iCount_ID]->str, pSegmentSlave, 4);
        str[iCount_ID]->str[4] = ' '; str[iCount_ID]->str[5] = pSegmentSlave[4];
        str[iCount_ID]->str[6] = ' '; str[iCount_ID]->str[7] = pSegmentSlave[5];
        str[iCount_ID]->str[8] = '\0';
        
        AlignId = pbsfThis->id;
        sprintf(AlignIdStr, "%d", AlignId);
        if(StringLen(AlignIdStr) > 8) {
           str[iCount_ID]->str[8] = ' '; str[iCount_ID]->str[9] = '\0';
           StringCat(str[iCount_ID]->str, AlignIdStr + 8);
        }           

        sip_cpy = MakeLocalId(str[iCount_ID]->str);
        iCount_ID++;

/*      oip = ObjectIdNew();
        oip->str = str; 

        sip = ValNodeNew(NULL);
        sip->choice = SEQID_LOCAL;
        sip->data.ptrvalue = oip; 

        sip_cpy = AsnIoMemCopy(sip, (AsnReadFunc) SeqIdAsnRead, (AsnWriteFunc) SeqIdAsnWrite);   */
       
        StringNCpy(PDBName, pSegmentSlave, 4);
        Chain = pSegmentSlave[4];

        AssignLocaIdToBiostruc(sip_cpy, PDBName, Chain, pbsSlaveHead); 
        AssignLocaIdToBioseq(sip_cpy, PDBName, Chain, sepHead);
        AssignLocaIdToSeqAlign(sip_cpy, PDBName, Chain, salpHead);     
     }   

     pbsfThis = pbsfThis->next;
     iCount1++;
     if(pbsSlaveHead) pbsSlaveHead = pbsSlaveHead->next;
     if(sepHead) sepHead = sepHead->next;
     if(salpHead) salpHead = salpHead->next;
  }

  return (pvnNcbi);

}
