/*   mmdbapi2.c
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
* File Name:  mmdbapi2.c
*
* Author:  Christopher Hogue, Hitomi Ohkawa
*
* Version Creation Date:   07/24/95
*
* $Revision: 6.3 $
*
* File Description:
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* 95/08/28 C. Hogue    Major change to new Biostruc ASN.1 v4.0,
*                      Addition of Features, Secondary Structure, Ensembles,
*                      PDB file generators, Kinemage file renderers moved
*                      to mmdbapi3.c
* 95/08/30 C. Hogue    Moved globals into mmdbapi2.c.
*
* $Log: mmdbapi2.c,v $
* Revision 6.3  1998/12/31 17:04:08  ywang
* work around reassign protein or NA as ion problem
*
 * Revision 6.2  1998/08/26  18:02:41  kans
 * fixed -v -fd warnings
 *
* Revision 6.1  1998/03/06 23:17:21  lewisg
* codewarrior fixes
*
* Revision 6.0  1997/08/25 18:11:03  madden
* Revision changed to 6.0
*
* Revision 5.5  1997/07/18 14:29:10  vakatov
* Provide return values (NULLs) in MolFromPDBChain() and SeqStringFromMol()
*
* Revision 5.4  1997/03/19 16:48:40  hogue
* Added fn SeqStringFromMol.
*
 * Revision 5.3  1997/01/27  21:41:10  hogue
 * Added function to find a chain based on PDB chain code (char)
 *
 * Revision 5.2  1996/08/19  21:08:10  vakatov
 * Removed redundant(and dangerous) function castings to "pNodeFunc"
 *
 * Revision 5.1  1996/07/15  21:28:00  hogue
 * Fixed glutamine labelling error.
 *
 * Revision 5.0  1996/05/28  14:02:09  ostell
 * Set to revision 5.0
 *
 * Revision 1.17  1996/05/22  20:43:16  hogue
 * Added new API functions GetFirstGraph, GetNextGraph, GetPrevGraph, GetMainAtom
 *
 * Revision 1.16  1996/05/06  21:23:14  hogue
 * added pointer tests to ToMSDParent, SetActiveModel, GetActiveModel
 *
 * Revision 1.15  1996/04/12  18:07:05  hogue
 * Changed flagging of atoms - AM_OCARBNYL now is anything branching off the backbone including N and H.  Added provisions for identifying D for Hydrogen in standard naming schemes.
 *
 * Revision 1.14  1996/03/29  19:34:57  hogue
 * Changed table of atomic radii, added H backbone flag, added ALD access in GetParent functions.
 *
 * Revision 1.13  1996/01/31  21:37:17  hogue
 * Removed extraneous traversing of structures - traversers limited to only one structure at a time.
 *
 * Revision 1.12  1995/11/15  19:16:37  hogue
 * TraverseSolids handles solids/densities on a per-model basis now.
 *
 * Revision 1.11  1995/10/10  12:27:08  hogue
 * Shut off unknown element error message
 *
 * Revision 1.10  1995/10/07  15:31:09  hogue
 * Selenium not an ion...
 *
 * Revision 1.9  1995/09/18  21:25:34  hogue
 * ProgMon() calls added.
 *
 * Revision 1.8  1995/09/14  15:56:00  kans
 * changed UInt1 to Uint1
 *
 * Revision 1.7  1995/09/14  14:22:30  hogue
 * Removed 254 and 255 in IsNotIon -> replaced with 0 (other) and INT1_MAX (unk)
 *
 * Revision 1.6  1995/09/13  21:50:32  hogue
 * Fixed Uint1 to array IsNotIon
 *
 * Revision 1.5  1995/09/05  19:10:33  hogue
 * Unknown element handling
 *
 * Revision 1.4  1995/08/30  19:36:44  kans
 * numerous changes
 *
 * Revision 1.3  1995/08/28  19:31:53  kans
 * version 4.1
 *
 * Revision 1.1  1995/08/02  17:07:15  kans
 * Initial revision
 *
*
* ==========================================================================
*/
/*****************************************************************************
*
*  mmdbapi2.c
*
*  UNDER CONSTRUCTION NOTICE.
*	SUBJECT TO RADICAL CHANGES...
*
*  programmer C.W.V. Hogue
*  created: 24 JULY 95
*  last mod: 28 AUGUST 95
*****************************************************************************/


#include <ncbi.h>
#include <asn.h>
#include <mmdbapi1.h>
#include <mmdbapi2.h>

CharPtr NCBIstdaaUC = "-ABCDEFGHIKLMNPQRSTVWXYZU*";
CharPtr NCBI4naUC = "-ACMGRSVTUWYHKDBN";
CharPtr NCBI4naLC = "-acmgrsvtuwyhkdbn";
CharPtr NCBIstdaaLC = "-abcdefghiklmnpqrstvwxyzu*";


/* these are acessed via AminoAcidName() */
static CharPtr NCBIstdaaMC3[MAX_NCBIstdaa] = {
	"---",
	"Ala",
	"Asx",
	"Cys",
	"Asp",
	"Glu",
	"Phe",
	"Gly",
	"His",
	"Ile",
	"Lys",
	"Leu",
	"Met",
	"Asn",
	"Pro",
	"Gln",
	"Arg",
	"Ser",
	"Thr",
	"Val",
	"Trp",
	"Xxx",
	"Tyr",
	"Glx",
	"Sec",
	"Ter" };

static  CharPtr NCBIstdaaUC3[MAX_NCBIstdaa] = {
	"---",
	"ALA",
	"ASX",
	"CYS",
	"ASP",
	"GLU",
	"PHE",
	"GLY",
	"HIS",
	"ILE",
	"LYS",
	"LEU",
	"MET",
	"ASN",
	"PRO",
	"GLN",
	"ARG",
	"SER",
	"THR",
	"VAL",
	"TRP",
	"XXX",
	"TYR",
	"GLX",
	"SEC",
	"TER" };


static  CharPtr NCBIstdaaLC3[MAX_NCBIstdaa] = {
	"---",
	"ala",
	"asx",
	"cys",
	"asp",
	"glu",
	"phe",
	"gly",
	"his",
	"ile",
	"lys",
	"leu",
	"met",
	"asn",
	"pro",
	"gln",
	"arg",
	"ser",
	"thr",
	"val",
	"trp",
	"xxx",
	"tyr",
	"glx",
	"sec",
	"ter" };


static CharPtr NCBIstdaaLCF[MAX_NCBIstdaa] = {
	"---",
	"alanine",
	"asx",
	"cysteine",
	"aspartate",
	"glutamate",
	"phenylalanine",
	"glycine",
	"histidine",
	"isoleucine",
	"lysine",
	"leucine",
	"methionine",
	"asparagine",
	"proline",
	"glutamine",
	"arginine",
	"serine",
	"threonine",
	"valine",
	"tryptophan",
	"xxx",
	"tyrosine",
	"glx",
	"selenocysteine",
	"termination" };


/* accessed via ElementName() */
static CharPtr IUPACElement[MAX_ELEMENTS] = {
 	" ",
	"H",
	"He",
	"Li",
	"Be",
	"B",
	"C",
	"N",
	"O",
	"F",
	"Ne",
	"Na",
	"Mg",
	"Al",
	"Si",
	"P",
	"S",
	"Cl",
	"Ar",
	"K",
	"Ca",
	"Sc",
	"Ti",
	"V",
	"Cr",
	"Mn",
	"Fe",
	"Co",
	"Ni",
	"Cu",
	"Zn",
	"Ga",
	"Ge",
	"As",
	"Se",
	"Br",
	"Kr",
	"Rb",
	"Sr",
	"Y",
	"Zr",
	"Nb",
	"Mo",
	"Tc",
	"Ru",
	"Rh",
	"Pd",
	"Ag",
	"Cd",
	"In",
	"Sn",
	"Sb",
	"Te",
	"I",
	"Xe",
	"Cs",
	"Ba",
	"La",
	"Ce",
	"Pr",
	"Nd",
	"Pm",
	"Sm",
	"Eu",
	"Gd",
	"Tb",
	"Dy",
	"Ho",
	"Er",
	"Tm",
	"Yb",
	"Lu",
	"Hf",
	"Ta",
	"W",
	"Re",
	"Os",
	"Ir",
	"Pt",
	"Au",
	"Hg",
	"Tl",
	"Pb",
	"Bi",
	"Po",
	"At",
	"Rn",
	"Fr",
	"Ra",
	"Ac",
	"Th",
	"Pa",
	"U",
	"Np",
	"Pu",
	"Am",
	"Cm",
	"Bk",
	"Cf",
	"Es",
	"Fm",
	"Md",
	"No",
	"Lr"};



/* this is from the Sargent-welch periodic table */
/* Bigger of atomic/covalent radius in Angstroms */
/* accessed via ElementSize() */
static FloatLo CommonSize[MAX_ELEMENTS] = {
 	0.0,
	0.79,
	0.93,
	2.05,
	1.40,
	1.17,
	1.30, /* fuued C */
	1.20, /* fuued N */
	1.10, /* fuued O */
	0.72,
	0.71,
	2.23,
	1.72,
	1.82,
	1.46,
	1.35, /* fuued P */
	1.30, /* fuued S */
	0.99,
	0.98,
	2.77,
	2.23,
	2.09,
	2.00,
	1.92,
	1.85,
	1.79,
	1.72,
	1.67,
	1.62,
	1.57,
	1.53,
	1.81,
	1.52,
	1.33,
	1.22,
	1.12,
	1.12,
	2.98,
	2.45,
	2.27,
	2.56,
	2.08,
	2.01,
	1.95,
	1.89,
	1.83,
	1.79,
	1.75,
	1.71,
	2.00,
	1.71,
	1.53,
	1.42,
	1.33,
	1.31,
	3.34,
	2.78,
	2.74,
	2.70,
	2.67,
	2.64,
	2.62,
	2.59,
	2.56,
	2.54,
	2.51,
	2.49,
	2.47,
	2.45,
	2.42,
	2.40,
	2.25,
	2.16,
	2.09,
	2.02,
	1.97,
	1.92,
	1.87,
	1.83,
	1.79,
	1.76,
	2.08,
	1.81,
	1.63,
	1.53,
	1.45,
	1.34,
	3.40, /* rest are est*/
	3.00,
	2.90,
	2.85,
	2.80,
	2.75,
	2.63,  /*rest are est */
	2.62,
	2.85,
	2.61,
	2.59,
	2.59,
	2.58,
	2.57,
	2.56,
	2.74,
	2.56};


static Uint1 IsNotIon[MAX_NOT_ION] = {0,1,2,5,6,7,8,10,14,15,16,18,34,36,54,86,INT1_MAX};


/* END OF GLOBAL VARIABLES */



Int2 LIBCALL TraverseSolids( DValNodePtr pdnModel, Int4 iModel, Int4 iIndex,
			  Pointer ptr, pNodeFunc pfnCallMe)
/*	 Walks the structure children from the current pointer.  Performs
	 (*pfnCallMe)((PFB) pxxxxThis,  Int4 iMode, Int4 iIndex, Pointer ptr) at MDD,  MSD
	 (as Flags Block). callbacks are void return types.

*/
{
    PFB pfbThis;
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    PVNMO pvnmoThis = NULL;
    PMOD pmodThis = NULL;
    PVNMD pvnmdThis = NULL;
    PMDD pmddThis = NULL;
    if (pdnModel)
      pfbThis = (PFB) pdnModel->data.ptrvalue;  /* get the data block */
    else return FALSE;
    switch ((int)pfbThis->bMe)  /* set the proper type */
      {
	  case AM_MSD:
	      pdnmsThis = (PDNMS) pdnModel;
	    break;
          default:
            return FALSE;
      }
    if (pdnmsThis)
     {
          pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
	  /* traverse Objects */
	  pvnmoThis = pmsdThis->pvnmoHead;
	  while (pvnmoThis)
	    {
 		if (pvnmoThis->choice == iModel)
  		  {	/* only those models stated */
			pmodThis = (PMOD) pvnmoThis->data.ptrvalue;
			(*pfnCallMe)((PFB) pmodThis, iModel, iIndex, ptr);
		  }
		pvnmoThis = pvnmoThis->next;
	    }

          /* traverse Densities */
	  pvnmdThis = pmsdThis->pvnmdHead;
	  while (pvnmdThis)
	    {
		if (pvnmdThis->choice == iModel)
		  {	/* only those models stated */
			pmddThis = (PMDD) pvnmdThis->data.ptrvalue;
			(*pfnCallMe)((PFB) pmddThis, iModel, iIndex, ptr);
		  }
		pvnmdThis = pvnmdThis->next;
	    }
     return TRUE;
    }
  return FALSE;
}



Int2 LIBCALL TraverseBonds( DValNodePtr pdnModel, Int4 iModel, Int4 iIndex,
			  Pointer ptr, pNodeFunc pfnCallMe)
/*	 Walks the structure children from the current pointer.  Performs
	 (*pfnCallMe)((PFB) pxxxxThis,  Int4 iMode, Int4 iIndex, Pointer ptr) at each PMBD
	 (as Flags Block). callbacks are void return types.

*/
{
    PFB pfbThis;
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    PDNMM pdnmmThis = NULL;
    PMMD pmmdThis = NULL;
    PDNMG pdnmgThis = NULL;
    PMGD pmgdThis = NULL;
    PVNMB pvnmbThis = NULL;
    PMBD pmbdThis = NULL;
    if (pdnModel)
      pfbThis = (PFB) pdnModel->data.ptrvalue;  /* get the data block */
    else return FALSE;
    switch ((int)pfbThis->bMe)  /* set the proper type */
      {
	  case AM_MSD:
	      pdnmsThis = (PDNMS) pdnModel;
	    break;
	  case AM_MMD:
	      pdnmmThis = (PDNMM) pdnModel;
	    break;
	  case AM_MGD:
	      pdnmgThis = (PDNMM) pdnModel;
	    break;
          default:
            return FALSE;
      }

    if (pdnmgThis)
       {
        while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      pvnmbThis = pmgdThis->pvnmbBHead;
		      while (pvnmbThis)
		        {   /*bond*/
			   pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
			   (*pfnCallMe)((PFB) pmbdThis,  iModel, iIndex, ptr); /* do the callback */
			   pvnmbThis = pvnmbThis->next;
			} /*bond*/

		      pdnmgThis = pdnmgThis->next;
		  }/*graph*/
        return TRUE;
       } /* if pdnmgThis */
    if (pdnmmThis)
       {
          while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
		pdnmgThis =  pmmdThis->pdnmgHead;
		while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      pvnmbThis = pmgdThis->pvnmbBHead;
		      while (pvnmbThis)
		        {   /*bond*/
			   pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
			   (*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr); /* do the callback */
			   pvnmbThis = pvnmbThis->next;
			} /*bond*/
		      pdnmgThis = pdnmgThis->next;
		  }/*graph*/
		/* traverse inter-residue bonds */
		pvnmbThis = pmmdThis->pvnmbIRBHead;
		while (pvnmbThis)
		  { /* inter-residue bonds */
		    pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
		    (*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr);  /* do the callback */
		    pvnmbThis = pvnmbThis->next;
		  }
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
        return TRUE;
       } /* if pdnmmThis */
    if (pdnmsThis)
     {
          pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
	  pdnmmThis = pmsdThis->pdnmmHead;
	  while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
	        pdnmgThis =  pmmdThis->pdnmgHead;
		while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      pvnmbThis = pmgdThis->pvnmbBHead;
		      while (pvnmbThis)
		        {   /*bond*/
			   pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
			   (*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr); /* do the callback */
			   pvnmbThis = pvnmbThis->next;
			} /*bond*/
		      pdnmgThis = pdnmgThis->next;
		  }/*graph*/
		/* traverse inter-residue bonds */
		pvnmbThis = pmmdThis->pvnmbIRBHead;
		while (pvnmbThis)
		  { /* inter-residue bonds */
		    pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
		    (*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr);  /* do the callback */
		    pvnmbThis = pvnmbThis->next;
		  }
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
	  /* traverse inter-molecule bonds */
	  pvnmbThis = pmsdThis->pvnmbIMBHead;
	  while (pvnmbThis)
	    { /* inter-molecule bonds */
		pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
		(*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr); /* do the callback */
		pvnmbThis = pvnmbThis->next;
	    }
     return TRUE;
    }
  return FALSE;
}




Int2 LIBCALL TraverseIBonds( DValNodePtr pdnModel, Int4 iModel, Int4 iIndex,
			  Pointer ptr, pNodeFunc pfnCallMe)
/*	 Walks the structure children from the current pointer.  Performs
	 (*pfnCallMe)((PFB) pxxxxThis,  Int4 iMode, Int4 iIndex, Pointer ptr)
	 at each PMBD Inter-residue or Inter-molecule (as Flags Block).
	 callbacks are void return types.

*/
{
    PFB pfbThis;
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    PDNMM pdnmmThis = NULL;
    PMMD pmmdThis = NULL;
    PVNMB pvnmbThis = NULL;
    PMBD pmbdThis = NULL;
    if (pdnModel)
      pfbThis = (PFB) pdnModel->data.ptrvalue;  /* get the data block */
    else return FALSE;
    switch ((int)pfbThis->bMe)  /* set the proper type */
      {
	  case AM_MSD:
	      pdnmsThis = (PDNMS) pdnModel;
	    break;
	  case AM_MMD:
	      pdnmmThis = (PDNMM) pdnModel;
	    break;
         default:
            return FALSE;
      }

    if (pdnmmThis)
       {
          while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
		/* traverse inter-residue bonds */
		pvnmbThis = pmmdThis->pvnmbIRBHead;
		while (pvnmbThis)
		  { /* inter-residue bonds */
		    pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
		    (*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr);  /* do the callback */
		    pvnmbThis = pvnmbThis->next;
		  }
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
        return TRUE;
       } /* if pdnmmThis */
    if (pdnmsThis)
     {
          pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
	  pdnmmThis = pmsdThis->pdnmmHead;
	  while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
		/* traverse inter-residue bonds */
		pvnmbThis = pmmdThis->pvnmbIRBHead;
		while (pvnmbThis)
		  { /* inter-residue bonds */
		    pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
		    (*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr);  /* do the callback */
		    pvnmbThis = pvnmbThis->next;
		  }
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
	  /* traverse inter-molecule bonds */
	  pvnmbThis = pmsdThis->pvnmbIMBHead;
	  while (pvnmbThis)
	    { /* inter-molecule bonds */
		pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
		(*pfnCallMe)((PFB) pmbdThis, iModel, iIndex, ptr); /* do the callback */
		pvnmbThis = pvnmbThis->next;
	    }
     return TRUE;
    }
  return FALSE;
}




Int2 LIBCALL TraverseAtoms( DValNodePtr pdnModel, Int4 iModel, Int4 iIndex,
			  Pointer ptr, pNodeFunc pfnCallMe)
/*	 Walks the structure children from the current pointer.  Performs
	 (*pfnCallMe)((PFB) pxxxxThis,  Int4 iMode, Int4 iIndex, Pointer ptr) at each PMAD
	 (as Flags Block). callbacks are void return types.

*/
{
    PFB pfbThis;
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    PDNMM pdnmmThis = NULL;
    PMMD pmmdThis = NULL;
    PDNMG pdnmgThis = NULL;
    PMGD pmgdThis = NULL;
    PVNMA pvnmaThis = NULL;
    PMAD pmadThis = NULL;
    if (pdnModel)
      pfbThis = (PFB) pdnModel->data.ptrvalue;  /* get the data block */
    else return FALSE;
    switch ((int)pfbThis->bMe)  /* set the proper type */
      {
	  case AM_MSD:
	      pdnmsThis = (PDNMS) pdnModel;
	    break;
	  case AM_MMD:
	      pdnmmThis = (PDNMM) pdnModel;
	    break;
	  case AM_MGD:
	      pdnmgThis = (PDNMM) pdnModel;
	    break;
          default:
            return FALSE;
      }

    if (pdnmgThis)
       {
        while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      pvnmaThis = pmgdThis->pvnmaAHead;
		      while (pvnmaThis)
		        {   /*atom*/
			    pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
			    (*pfnCallMe)((PFB) pmadThis, iModel, iIndex, ptr); /* do the callback */

			    pvnmaThis = pvnmaThis->next;
			}   /*atom*/
		       pdnmgThis = pdnmgThis->next;
		  }/*graph*/
        return TRUE;
       } /* if pdnmgThis */
    if (pdnmmThis)
       {
          while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
	        pdnmgThis =  pmmdThis->pdnmgHead;
		while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      pvnmaThis = pmgdThis->pvnmaAHead;
		      while (pvnmaThis)
		        {   /*atom*/
			    pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
			    (*pfnCallMe)((PFB) pmadThis, iModel, iIndex, ptr); /* do the callback */

			    pvnmaThis = pvnmaThis->next;
			}   /*atom*/
		      pdnmgThis = pdnmgThis->next;
		  }/*graph*/
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
        return TRUE;
       } /* if pdnmmThis */
    if (pdnmsThis)
     {
          pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
          pdnmmThis = pmsdThis->pdnmmHead;
	  while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
	        pdnmgThis =  pmmdThis->pdnmgHead;
		while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      pvnmaThis = pmgdThis->pvnmaAHead;
		      while (pvnmaThis)
		        {   /*atom*/
			    pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
			    (*pfnCallMe)((PFB) pmadThis, iModel, iIndex, ptr); /* do the callback */

			    pvnmaThis = pvnmaThis->next;
			}   /*atom*/

		      pdnmgThis = pdnmgThis->next;
		  }/*graph*/
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
     return TRUE;
    }
  return FALSE;
}




Int2 LIBCALL TraverseGraphs( DValNodePtr pdnModel, Int4 iModel, Int4 iIndex,
			  Pointer ptr, pNodeFunc pfnCallMe)
/*	 Walks the structure children from the current pointer.  Performs
	 (*pfnCallMe)((PFB) pxxxxThis,  Int4 iMode, Int4 iIndex, Pointer ptr) at PMGD
	 (as Flags Block). callbacks are void return types.

*/
{
    PFB pfbThis;
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    PDNMM pdnmmThis = NULL;
    PMMD pmmdThis = NULL;
    PDNMG pdnmgThis = NULL;
    PMGD pmgdThis = NULL;
    if (pdnModel)
      pfbThis = (PFB) pdnModel->data.ptrvalue;  /* get the data block */
    else return FALSE;
    switch ((int)pfbThis->bMe)  /* set the proper type */
      {
	  case AM_MSD:
	      pdnmsThis = (PDNMS) pdnModel;
	    break;
	  case AM_MMD:
	      pdnmmThis = (PDNMM) pdnModel;
	    break;
	  case AM_MGD:
	      pdnmgThis = (PDNMM) pdnModel;
	    break;
          default:
            return FALSE;
      }

    if (pdnmgThis)
       {
        while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      (*pfnCallMe)((PFB) pmgdThis, iModel, iIndex, ptr); /* do the callback */
		      pdnmgThis = pdnmgThis->next;
		  }/*graph*/
        return TRUE;
       } /* if pdnmgThis */
    if (pdnmmThis)
       {
          while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
	        pdnmgThis =  pmmdThis->pdnmgHead;
		while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      (*pfnCallMe)((PFB) pmgdThis, iModel, iIndex, ptr); /* do the callback */
		       pdnmgThis = pdnmgThis->next;
		  }/*graph*/
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
        return TRUE;
       } /* if pdnmmThis */
    if (pdnmsThis)
     {
          pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
	  pdnmmThis = pmsdThis->pdnmmHead;
	  while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
	        pdnmgThis =  pmmdThis->pdnmgHead;
		while(pdnmgThis)
		  {   /*graph*/
		      pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		      (*pfnCallMe)((PFB) pmgdThis, iModel, iIndex, ptr); /* do the callback */
		      pdnmgThis = pdnmgThis->next;
		  }/*graph*/
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
     return TRUE;
    }
  return FALSE;
}




Int2 LIBCALL TraverseMolecules( DValNodePtr pdnModel, Int4 iModel, Int4 iIndex,
			  Pointer ptr, pNodeFunc pfnCallMe)
/*	 Walks the structure children from the current pointer.  Performs
	 (*pfnCallMe)((PFB) pxxxxThis,  Int4 iMode, Int4 iIndex, Pointer ptr) at each PMMD
	 (as Flags Block). callbacks are void return types.

*/
{
    PFB pfbThis;
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    PDNMM pdnmmThis = NULL;
    PMMD pmmdThis = NULL;

    if (pdnModel)
      pfbThis = (PFB) pdnModel->data.ptrvalue;  /* get the data block */
    else return FALSE;
    switch ((int)pfbThis->bMe)  /* set the proper type */
      {
	  case AM_MSD:
	      pdnmsThis = (PDNMS) pdnModel;
	    break;
	  case AM_MMD:
	      pdnmmThis = (PDNMM) pdnModel;
	    break;
          default:
            return FALSE;
      }
    if (pdnmmThis)
       {
          while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
		(*pfnCallMe)((PFB) pmmdThis, iModel, iIndex, ptr); /* do the callback */
		pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
        return TRUE;
       } /* if pdnmmThis */
    if (pdnmsThis)
     {
          pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
	  pdnmmThis = pmsdThis->pdnmmHead;
	  while(pdnmmThis)
	    {   /*molecule*/
		pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
		(*pfnCallMe)((PFB) pmmdThis, iModel, iIndex, ptr); /* do the callback */
	      pdnmmThis = pdnmmThis->next;
	    }   /*molecule*/
     return TRUE;
    }
  return FALSE;
}



Int2 LIBCALL TraverseModels(DValNodePtr pdnModel, Int2 iTraverse,
     Int4 iIndex, Pointer ptr,  pNodeFunc pfnCallMe)
{
  /* this routine iterates through all the models in a MSD */
  /* it sets the iActiveModel setting in the parent MSD    */
  /* then calls TraverseAll - it is up to you to look at   */
  /* iActiveModel with GetActiveModel(pfbThis) */
  /* you are allowed to call this from the level of */
  /* structure, molecule or graph */
  PMSD pmsdThis = NULL;
  PDNML pdnmlThis = NULL;
  Int2 iHoldActive = 1;
  Int2 iTest = FALSE;
  pmsdThis = ToMSDParent(PFBFromDN(pdnModel));
  if (!pmsdThis) return 0;
  pdnmlThis = pmsdThis->pdnmlModels; /* the linked-list of models */
  iHoldActive = pmsdThis->iActiveModel;
  if (pdnmlThis)
  while (pdnmlThis)
    {
      pmsdThis->iActiveModel = (Int2) pdnmlThis->choice;
      switch (iTraverse)
        {
	    case TRAVERSE_SOLID:
	       iTest = TraverseSolids( pdnModel, (Int4) pdnmlThis->choice, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_BOND:
	       iTest = TraverseBonds( pdnModel, (Int4) pdnmlThis->choice, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_IBOND:
	       iTest = TraverseIBonds( pdnModel, (Int4) pdnmlThis->choice, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_ATOM:
	       iTest = TraverseAtoms( pdnModel, (Int4) pdnmlThis->choice, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_GRAPH:
	       iTest = TraverseGraphs( pdnModel, (Int4) pdnmlThis->choice, iIndex, ptr, pfnCallMe);
               break;
	    case TRAVERSE_MOLECULE:
	       iTest = TraverseMolecules( pdnModel, (Int4) pdnmlThis->choice, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_ALL:
	    default:
	       iTest = TraverseAll( pdnModel, (Int4) pdnmlThis->choice, iIndex, ptr, pfnCallMe);

	}
      if (iTest == FALSE) break;
      pdnmlThis = pdnmlThis->next;
    }
  pmsdThis->iActiveModel = iHoldActive;
  return iTest;
}


Int2 LIBCALL TraverseOneModel(DValNodePtr pdnModel, Int2 iTraverse, Int2 iModel,
			Int4 iIndex, Pointer ptr, pNodeFunc pfnCallMe)
{
  Int2 iHoldActive = 1;
  Int2 iTest = FALSE;

  iHoldActive = SetActiveModel(PFBFromDN(pdnModel), iModel);
     switch (iTraverse)
        {
	    case TRAVERSE_SOLID:
	       iTest = TraverseSolids( pdnModel, (Int4) iModel, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_BOND:
	       iTest = TraverseBonds( pdnModel, (Int4) iModel, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_IBOND:
	       iTest = TraverseIBonds( pdnModel, (Int4) iModel, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_ATOM:
	       iTest = TraverseAtoms( pdnModel, (Int4) iModel, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_GRAPH:
	       iTest = TraverseGraphs( pdnModel, (Int4) iModel, iIndex, ptr, pfnCallMe);
               break;
	    case TRAVERSE_MOLECULE:
	       iTest = TraverseMolecules( pdnModel, (Int4) iModel, iIndex, ptr, pfnCallMe);
	       break;
	    case TRAVERSE_ALL:
	    default:
	       iTest = TraverseAll( pdnModel, (Int4) iModel, iIndex, ptr, pfnCallMe);

	}
  iModel = SetActiveModel(PFBFromDN(pdnModel), iHoldActive);
  return iTest;
}





Int2 LIBCALL IndexFromNode(PFB pfbNode)
{
    PMSD  pmsdThis;
    PMMD  pmmdThis;
    PMGD  pmgdThis;
    PMAD  pmadThis;
    if (!pfbNode) return 0;
    switch ((int)pfbNode->bMe)
      {
	  case AM_MSD:
		pmsdThis = (PMSD) pfbNode;
		return pmsdThis->pdnmsLink->choice;
	  case AM_MMD:
		pmmdThis = (PMMD) pfbNode;
		return pmmdThis->pdnmmLink->choice;
	  case AM_MGD:
		pmgdThis = (PMGD) pfbNode;
		return pmgdThis->pdnmgLink->choice;
	  case AM_MAD:
	        pmadThis = (PMAD) pfbNode;
		return pmadThis->iIndex;
	  case AM_MBD:
	  case AM_MOD:
          case AM_MDD:
	  case AM_UNK:
	     return 0;
      }
    return 0;
}

PMSD LIBCALL ToMSDParent(PFB pfbThis)
 {
   PMSD pmsdThis = NULL;

    if (! pfbThis) return NULL;
     if (!pfbThis->pfbParent)
       {
         if (pfbThis->bMe & (Byte) AM_MSD) return (PMSD) pfbThis;
       }
     else
       {
        pmsdThis = ToMSDParent((PFB) pfbThis->pfbParent);
	return pmsdThis;
       }
    return NULL;
 }

CharPtr LIBCALL ParentStrucName(PFB pfbThis)
 {
   PMSD pmsdThis = NULL;
   pmsdThis = ToMSDParent(pfbThis);
   if (pmsdThis) return pmsdThis->pcPDBName;
    else return NULL;
 }

PMGD LIBCALL GetParentGraph(PFB pfbThis)
 {  /* use for atoms and bonds */
     PFB pfbHold = NULL;

     pfbHold = pfbThis->pfbParent;
     if (pfbThis->bMe == (Byte) AM_ALD) pfbHold = pfbHold->pfbParent;
     if ((pfbThis->bMe & (Byte) AM_MAD) || (pfbThis->bMe & (Byte) AM_MBD)
         || (pfbThis->bMe & (Byte) AM_ALD))
       {
         if (pfbHold->bMe & (Byte) AM_MGD)
	    return (PMGD) pfbHold;
	 else return NULL;  /* IRBonds and IMBonds have no parent ! */
       }
     else return NULL;
 }

CharPtr LIBCALL ParentGraphPDBNo(PFB pfbThis)
 {
     PMGD pmgdThis = NULL;
     pmgdThis = GetParentGraph(pfbThis);
     if (pmgdThis)
        return pmgdThis->pcGraphNum;
     else return NULL;
 }

CharPtr LIBCALL ParentGraphPDBName(PFB pfbThis)
 {
     PMGD pmgdThis = NULL;
     pmgdThis = GetParentGraph(pfbThis);
     if (pmgdThis)
        return pmgdThis->pcGraphName;
     else return NULL;
 }

CharPtr LIBCALL ParentGraphIUPAC1(PFB pfbThis)
 {  /* returns X when not one of 20 amino acids or 4 NT bases */
     PMGD pmgdThis = NULL;
     pmgdThis = GetParentGraph(pfbThis);
     if (pmgdThis)
        return pmgdThis->pcIUPAC;
     else return NULL;
 }


Int2 LIBCALL ParentGraphNum(PFB pfbThis)
 {
   PMGD pmgdThis = NULL;
     pmgdThis = GetParentGraph(pfbThis);
   if (pmgdThis)
      {
        return pmgdThis->pdnmgLink->choice;
      }
   return 0;
 }

PMMD LIBCALL GetParentMol(PFB pfbThis)
 { /* use for bonds,  atom,  graphs */
   PFB pfbHold = NULL;
   if (!pfbThis) return NULL;
   pfbHold = pfbThis->pfbParent;
   switch ((int)pfbThis->bMe)
    {
      case AM_MBD:
          switch ((int)pfbHold->bMe)
	     {
	       case AM_MGD:  /* graph bond */
	           pfbHold = pfbHold->pfbParent;
		   if (pfbHold->bMe & (Byte) AM_MMD) return (PMMD) pfbHold;
		   else return NULL;
	       case AM_MMD:  /* inter-res Bond */
	           return (PMMD) pfbHold;
	       case AM_MSD:  /* inter-mol bond */
	       default:
		    return NULL; /* no parent */
	     }
       case AM_ALD:
           pfbHold = pfbHold->pfbParent;
       case AM_MAD:
	   if (pfbHold->bMe & (Byte) AM_MGD) pfbHold = pfbHold->pfbParent;
	   else return NULL;
	   if (pfbHold->bMe & (Byte) AM_MMD) return (PMMD) pfbHold;
	   else return NULL;
       case AM_MGD:
	   if (pfbHold->bMe & (Byte) AM_MMD) return (PMMD) pfbHold;
	   else return NULL;
       default:
	   pfbHold = NULL;
    }
    return (PMMD) pfbHold;
 }




CharPtr LIBCALL ParentMolName(PFB pfbThis)
 {
     PMMD pmmdThis;
     pmmdThis = GetParentMol(pfbThis);
     if (pmmdThis) return pmmdThis->pcMolName;
     else return NULL;
 }


Int2 LIBCALL ParentMolNum(PFB pfbThis)
 {
   PMMD pmmdThis = NULL;
   pmmdThis = GetParentMol(pfbThis);
   if (pmmdThis)
      {
        return pmmdThis->pdnmmLink->choice;
      }
   return 0;
 }

PFB LIBCALL PFBFromVN(ValNodePtr pvnModel)
 {
     PFB pfbThis = NULL;
     pfbThis =   (PFB) pvnModel->data.ptrvalue;
     return pfbThis;
 }

PFB LIBCALL PFBFromDN(DValNodePtr pdnModel)
 {
     PFB pfbThis = NULL;
     pfbThis =  (PFB) pdnModel->data.ptrvalue;
     return pfbThis;
 }

DValNodePtr DNFromPFB(PFB pfbThis)
{
    PMSD pmsdThis;
    PMMD pmmdThis;
    PMGD pmgdThis;
    DValNodePtr pdnThis = NULL;
    if (pfbThis)
    switch ((int)pfbThis->bMe)
      {
	case AM_MSD:
	    pmsdThis = (PMSD) pfbThis;
	    pdnThis = (DValNodePtr) pmsdThis->pdnmsLink;
	    return pdnThis;
	case AM_MMD:
	    pmmdThis = (PMMD) pfbThis;
	    pdnThis = (DValNodePtr) pmmdThis->pdnmmLink;
	    return pdnThis;
	case AM_MGD:
	    pmgdThis = (PMGD) pfbThis;
	    pdnThis = (DValNodePtr) pmgdThis->pdnmgLink;
	    return pdnThis;
	default:
	  return NULL;
      }
   return NULL;
}


Int2 LIBCALL SetActiveModel(PFB pfbThis, Int2 iModel)
 {
   PMSD pmsdThis;
   Int2 iOldModel;
   if (!pfbThis) return 0;
   pmsdThis = ToMSDParent(pfbThis);
   iOldModel = pmsdThis->iActiveModel;
   pmsdThis->iActiveModel = iModel;
   return iOldModel;
 }


Int2 LIBCALL GetActiveModel(PFB pfbThis)
 {
     PMSD pmsdThis = NULL;
     if (!pfbThis) return 0;
     pmsdThis = (PMSD) ToMSDParent(pfbThis);
     return pmsdThis->iActiveModel;
 }


PALD LIBCALL GetAtomLocs(PMAD pmadThis,  Int2 iModel)
 {
     PVNAL pvnalThis = NULL;
     if (pmadThis)
     if (iModel == 0)
	{
	   if (pmadThis->pvnalLocate)
	     return (PALD) pmadThis->pvnalLocate->data.ptrvalue;
	}
     else
	{
	    pvnalThis = ValNodeFindNext(pmadThis->pvnalLocate, NULL, iModel);
	     if (pvnalThis) return (PALD) pvnalThis->data.ptrvalue;
	}
    return NULL;
 }




/* can call: paldThis = GetLocation(GetAtomLocs(pmadThis, GetActiveModel(pmadThis)), 0, ' ') */
/* need a GetActiveCoordSet(iModel), GetEnsembleString(iModel, iCoordset ) */




PALD LIBCALL GetLocation(PALD paldList,  Int2 iCoord,  Char cAlt)
{
    PALD paldThis = NULL;
    if (iCoord) ; /* later eh */
    paldThis = paldList;
    while (paldThis)
       {
	    if (paldThis->cAltConf == cAlt)
	         {
		     return paldThis;
		 }
	    paldThis = paldThis->next;
       }
    return NULL;
}





/* not for external-use only for AssignAtomLocId */
static void LIBCALLBACK UniqueLocId(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr)
{
   PMAD pmadThis = NULL;
   PALD paldThis = NULL;
   Int4Ptr piIndex;
   piIndex = (Int4Ptr) ptr;
   if (pfbThis->bMe == (Byte) AM_MAD)
     {  /* iIndex isn't used */
	pmadThis = (PMAD) pfbThis;
	paldThis = GetAtomLocs(pmadThis, iModel);
 	while (paldThis)
	   {  /* dosen't care about coord-set id's, alt-confs or ensembles */
	     *piIndex = (Int4) (*piIndex + 1);
	     paldThis->iUniqueId = (Int4) (*piIndex);
	     paldThis = paldThis->next;
	   }
     }
}


Int4 LIBCALL AssignAtomLocId(PDNMS pdnmsThis)
{
    Int2 iTest = 0;
    Int4 iReservedIndex = 0;
    ProgMon("Assigning Location ID's");
    iTest = TraverseModels(pdnmsThis, TRAVERSE_ATOM,
              0, (Pointer) &iReservedIndex, UniqueLocId );
    return iReservedIndex;  /* the number of locs assigned - also the number of xyz coords */
}


CharPtr LIBCALL ElementName(Int1 iAtomicNo)
{
	if ((iAtomicNo < MAX_ELEMENTS) && (iAtomicNo > 0))
	  return IUPACElement[iAtomicNo];
	else return NULL;
}

FloatLo LIBCALL ElementSize(Int1 iAtomicNo)
{
	if ((iAtomicNo < MAX_ELEMENTS) && (iAtomicNo > 0))
	  return  CommonSize[iAtomicNo];
	else return (FloatLo) 0.0;
}


CharPtr LIBCALL AminoAcidName(Char cThis, Int2 iRetCode)
{
	CharPtr pcTemp = NULL;
	Int2 iCount = 0;

	pcTemp = NCBIstdaaUC;
	while (pcTemp[iCount])
   	  {
		if (pcTemp[iCount] == cThis)
		  switch(iRetCode)
		   {
			case USE_UPPERCASE:
			  return NCBIstdaaUC3[iCount];
			case USE_LOWERCASE:
			  return NCBIstdaaLC3[iCount];
		 	case USE_FULLNAME:
			  return NCBIstdaaLCF[iCount];
			case USE_MIXCASE:
			default:
			  return NCBIstdaaMC3[iCount];
		   }
		iCount++;
	  }
	return NULL;
}


CharPtr LIBCALL AminoAcidNameFromIdx(Int2 iIndex, Int2 iRetCode)
{
	 if (iIndex < MAX_NCBIstdaa)
	    switch(iRetCode)
		   {
			case USE_UPPERCASE:
			  return NCBIstdaaUC3[iIndex];
			case USE_LOWERCASE:
			  return NCBIstdaaLC3[iIndex];
		 	case USE_FULLNAME:
			  return NCBIstdaaLCF[iIndex];
			case USE_MIXCASE:
			default:
			  return NCBIstdaaMC3[iIndex];
		   }
	return NULL;
}

void LIBCALLBACK SetBackBone(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr)
{
   PMAD pmadThis = NULL;
   PMMD pmmdThis = NULL;
   if (pfbThis->bMe == (Byte) AM_MAD)
     {  /* iIndex isn't used */
	pmadThis = (PMAD) pfbThis;
	if (pmadThis->pcAName)
	 {
	  pmmdThis = GetParentMol((PFB)pmadThis);
	  if ((pmmdThis->bWhat & (Byte) AM_DNA) || (pmmdThis->bWhat & (Byte) AM_RNA))
	    {
	      if (!StringICmp(pmadThis->pcAName, " P  "))
	         {
		    pmadThis->bWhat = (Byte)(pmadThis->bWhat |(Byte) AM_PALPHA | (Byte) AM_BACKBONE);
		    return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " O2P"))
	         {
		    pmadThis->bWhat = (Byte)(pmadThis->bWhat |(Byte) AM_OCARBNYL| (Byte) AM_BACKBONE);
		    return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " O1P"))
	         {
		    pmadThis->bWhat = (Byte)(pmadThis->bWhat |(Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		    return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " C1*"))
	         {
		    pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_C1RIBOSE );
		    return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " C4*"))
	         {
		    pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_C4RIBOSE | (Byte) AM_BACKBONE);
		    return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " C5*"))
	         {
		    pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_BACKBONE);
		    return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " O5*"))
	         {
		    pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_BACKBONE);
		    return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " C3*"))
	         {
		    pmadThis->bWhat = (Byte) (pmadThis->bWhat | (Byte) AM_BACKBONE);
		    return;
		 }
	       if (!StringICmp(pmadThis->pcAName, " O3*"))
	         {
		    pmadThis->bWhat = (Byte) (pmadThis->bWhat | (Byte) AM_BACKBONE);
		    return;
		 }
	       if (!StringICmp(pmadThis->pcAName, " N9 "))
	         {
		    pmadThis->bWhat = (Byte) (pmadThis->bWhat | (Byte) AM_NBETA);
		    return;
		 }
  /*     code that must test parent graph for Pyr or Pur
		if (!StringICmp(pmadThis->pcAName, " N1 "))
	         {    if A or G don't set this
		    pmadThis->bWhat = (Byte) (pmadThis->bWhat | (Byte) AM_NBETA);
		    return;
		 }
		*/
	      return;
	    }
	  if (pmmdThis->bWhat & (Byte) AM_PROT)
	     {
	       if (!StringICmp(pmadThis->pcAName," CA "))
		 {
		     pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_CALPHA | (Byte) AM_BACKBONE);
		     return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " C  "))
	         {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte)AM_BACKBONE);
		   return;
		 }
	       if (!StringICmp(pmadThis->pcAName, " N  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte)AM_BACKBONE);
		   return;
		 }
	      if (!StringICmp(pmadThis->pcAName, " H  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
              if (!StringICmp(pmadThis->pcAName, " D  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
               if (!StringICmp(pmadThis->pcAName, "1H  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
              if (!StringICmp(pmadThis->pcAName, "2H  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
 	       if (!StringICmp(pmadThis->pcAName, "3H  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
	       if (!StringICmp(pmadThis->pcAName, "1D  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
              if (!StringICmp(pmadThis->pcAName, "2D  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
 	       if (!StringICmp(pmadThis->pcAName, "3D  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
 	      if (!StringICmp(pmadThis->pcAName, " O  "))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
   	      if (!StringICmp(pmadThis->pcAName, " OXT"))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
              if (!StringICmp(pmadThis->pcAName, " HXT"))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
              if (!StringICmp(pmadThis->pcAName, " DXT"))
		 {
	           pmadThis->bWhat = (Byte)(pmadThis->bWhat | (Byte) AM_OCARBNYL | (Byte) AM_BACKBONE);
		   return;
		 }
	      if (!StringICmp(pmadThis->pcAName," CB "))
		 {
		     pmadThis->bWhat = (Byte) (pmadThis->bWhat | (Byte) AM_CBETA);
		     return;
		 }  /* dosen't do GLY */
	       return;
	     } /* if AM_PROT */
	 } /* if pcAName */
     }  /* if AM_MAD */
   return;
}


void LIBCALL AssignBackBone(PDNMS pdnmsThis)
{
    Int2 iTest = 0;
    ProgMon("Assigning Backbone");
    iTest =  TraverseAtoms(pdnmsThis,  0, 0, NULL,  SetBackBone);
    return;
}

FloatLo LIBCALL AtomDistanceSq(PMAD pmadFrom, PMAD pmadTo, Int2 iModel, Char  cAlt)
{
  PALD paldFrom = NULL;
  PALD paldTo = NULL;
  FloatLo flDist = 0.0;
  FloatLo flX, flY, flZ;


  if (!(pmadFrom)) return 0.0;
  if (!(pmadTo)) return 0.0;
  paldFrom = GetLocation(GetAtomLocs(pmadFrom, iModel), 0, cAlt);
  paldTo = GetLocation(GetAtomLocs(pmadTo, iModel), 0, cAlt);
  if ((!paldFrom) || (!paldTo)) return -1.0;
  flX = paldFrom->pflvData[0] - paldTo->pflvData[0];
  flY = paldFrom->pflvData[1] - paldTo->pflvData[1];
  flZ = paldFrom->pflvData[2] - paldTo->pflvData[2];
  flDist = flX * flX + flY * flY + flZ * flZ;
  return flDist;
}


void LIBCALLBACK LinkAlpha(PFB pfbThis,  Int4 iModel,  Int4 iIndex,  Pointer ptr)
{  /* only used while walking within a single molecule */
   /* iModel is NOT passed through */
    PMAD pmadThis = NULL;
    PMAD pmadHold = NULL;
    PPFB ppfbHold = NULL;
    PVNMB pvnmbThis = NULL;
    PMBD pmbdThis = NULL;
    PMGD pmgdThis = NULL;
    FloatLo flDistSq = 0.0;
    ValNodePtr pvnThis = NULL;
    if (pfbThis->bMe == (Byte) AM_MAD)
    if ((pfbThis->bWhat & (Byte) AM_CALPHA)  ||
	 (pfbThis->bWhat & (Byte) AM_PALPHA))
      {  /* do the link to the last Alpha pointed to by ptr */
  	ppfbHold = (PPFB) ptr;
	pmadHold = (PMAD) *ppfbHold;
	pmadThis = (PMAD) pfbThis;
        /* set the hold */
	*ppfbHold = pfbThis;
	if (pmadHold)
	  {  /* we have a bond */
	     pmgdThis = GetParentGraph(pfbThis);
	     /* check the distance for the outrageous */
	     flDistSq = AtomDistanceSq(pmadHold, pmadThis, GetActiveModel(pfbThis), ' ');
	     if ((pmgdThis->bWhat & RES_AA) && ( flDistSq > VIRTUAL_BOND_AA_MAX_SQ))
                 return;
       	     if (((pmgdThis->bWhat & RES_RNA) || (pmgdThis->bWhat & RES_DNA)) &&
	           ( flDistSq > VIRTUAL_BOND_NA_MAX_SQ))
                 return;

	     pvnmbThis = NewVNMB(&pmgdThis->pvnmbBHead,(Int2) VIRTUAL_BOND_CHOICE); /* 7 = Other */
	     if (!pvnmbThis) return;  /* callback fails silently */
	     pmgdThis->iBondCount += 1;
	     pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
	     pmbdThis->pfbParent = (PFB) pmgdThis;
	     pmbdThis->bWhat = (Byte)(pmbdThis->bWhat | (Byte) BOND_VIRTUAL);
	     pmbdThis->pmadFrom = pmadThis;
	     pvnThis = ValNodeAddPointer(&pmadThis->pvnBonds,
					 (Int2) VIRTUAL_BOND_CHOICE,
					 (Nlm_VoidPtr) pmbdThis);
	     if (!pvnThis) return;
	     pmbdThis->pmadTo = pmadHold;
	     pvnThis = ValNodeAddPointer(&pmadHold->pvnBonds,
					(Int2) VIRTUAL_BOND_CHOICE,
					(Nlm_VoidPtr) pmbdThis);
	     if (!pvnThis) return;
	  }
     }

   return;
}


void LIBCALLBACK SetVirtualBB(PFB pfbThis,  Int4 iModel,  Int4 iIndex,  Pointer ptr)
{
    PDNMG pdnmgThis = NULL;
    PMMD pmmdThis = NULL;
    PFB pfbHold = NULL;

    if (pfbThis->bMe == (Byte) AM_MMD)
      {
	  pmmdThis = (PMMD) pfbThis;
	  if ((pmmdThis->bWhat == (Byte) AM_PROT) ||
	     (pmmdThis->bWhat == (Byte) AM_DNA) ||
	      (pmmdThis->bWhat == (Byte) AM_RNA))
	    {
		pdnmgThis = pmmdThis->pdnmgHead;
		TraverseAtoms(pdnmgThis, iModel, iIndex, (Pointer) &pfbHold,
                   LinkAlpha);
	    }
      }
   return;
}

void LIBCALL AssignVirtual(PDNMS pdnmsThis)
{
    Int2 iTest = 0;
    ProgMon("Assigning Virtual Bonds");
    iTest = TraverseMolecules(pdnmsThis,  0,  0,  NULL,  SetVirtualBB);
   /* iTest = TraverseIBonds(pdnmsThis,  0, 0, NULL, SetHetVirt); */

}

void LIBCALLBACK SetIons(PFB pfbThis,  Int4 iModel,  Int4 iIndex,  Pointer ptr)
{
    PMAD pmadThis = NULL;
    PMMD pmmdThis = NULL;
    PMSD pmsdThis = NULL;
    Int2 iCount;
    Boolean bIon = FALSE;
    if (pfbThis->bMe == (Byte) AM_MAD)
      {
	  pmadThis = (PMAD) pfbThis;
	  pmmdThis = GetParentMol(pfbThis);
	  pmsdThis = ToMSDParent(pfbThis);
	  if ((Int1) pmadThis->pvnmaLink->choice > MAX_ELEMENTS)
	     {  /* report on unknown elements */
	/*	 MsgAlert(KEY_NONE,SEV_ERROR,"MMDB Warning", "Unknown Element\nNotify NCBI of Accession No: %s",pmsdThis->pcPDBName);
	*/	 return;
	     }
	  /* catch-C,O,N,H - will return faster than the other way */
	  /* for the large majority of biomolecular atoms */
	  if ((Int1) pmadThis->pvnmaLink->choice == 1) return;
	  if ((Int1) pmadThis->pvnmaLink->choice == 6) return;
	  if ((Int1) pmadThis->pvnmaLink->choice == 8) return;
	  if ((Int1) pmadThis->pvnmaLink->choice == 7) return;
	  bIon = TRUE; /* assume that it is an ion  */
	  for (iCount = 0; iCount < MAX_NOT_ION; iCount++)
	     {
	       if ((Uint1) pmadThis->pvnmaLink->choice == IsNotIon[iCount])
	         {  /* found it is not an ion */
		   bIon = FALSE;
		   break;
		 }
	     }
	  if (bIon) {
         if(pmmdThis->bWhat != AM_PROT && pmmdThis->bWhat != AM_DNA && pmmdThis->bWhat != AM_RNA){
		    pmmdThis->bWhat = (Byte) (AM_ION);
		    pmadThis->bUpdate = (Byte) (pmadThis->bUpdate | AM_ION);
		    pmsdThis->bWhat = (Byte) (pmsdThis->bWhat | AM_ION);
	     }
      }
      }
   return;
}





void LIBCALL AssignIons(PDNMS pdnmsThis)
{
    Int2 iTest = 0;
    ProgMon("Assigning Ions");
    iTest = TraverseAtoms(pdnmsThis,  0,  0,  NULL,  SetIons);
}



PDNMG LIBCALL GetFirstGraph(PMMD pmmdThis)
{
  if (!pmmdThis) return NULL;
  return pmmdThis->pdnmgHead;
}

PDNMG LIBCALL GetLastGraph(PMMD pmmdThis)
{
  PDNMG pdnmgThis = NULL;
  if (!pmmdThis) return NULL;
  pdnmgThis = pmmdThis->pdnmgHead;
  while (pdnmgThis)
   {
      if (pdnmgThis->next == NULL) return pdnmgThis;
      pdnmgThis = pdnmgThis->next;
   }
 return NULL;
}


PDNMG LIBCALL GetNextGraph(PDNMG pdnmgThis)
{
  if (!pdnmgThis) return NULL;
  return pdnmgThis->next;
}

PDNMG LIBCALL GetPrevGraph(PDNMG pdnmgThis)
{
  if (!pdnmgThis) return NULL;
  return pdnmgThis->last;
}

PMAD LIBCALL GetMainAtom(PDNMG pdnmgThis)
{
  PMGD pmgdThis = NULL;
  PVNMA pvnmaThis = NULL;
  PMAD pmadThis = NULL;
  Byte MainAtom = 0;

  if (!pdnmgThis) return NULL;
  pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
  if ((pmgdThis->bWhat & (Byte) RES_RNA) || (pmgdThis->bWhat & (Byte) RES_DNA))
     MainAtom = AM_PALPHA;
  if (pmgdThis->bWhat & (Byte) RES_AA) MainAtom = AM_CALPHA;
  pvnmaThis = pmgdThis->pvnmaAHead;
  while (pvnmaThis)
   {
     pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
     if (pmadThis->bWhat & MainAtom)
       return pmadThis;
     pvnmaThis = pvnmaThis->next;
   }
 return NULL;
}

void LIBCALLBACK FindMolByChar(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PMMD pmmdMol;
  Char c;
  PDNMM *ppdnmmFound;
  
    ppdnmmFound = (PDNMM *)ptr;
    if (pfbThis->bMe == AM_MMD)
      {
          /* cast to the correct PMMD pointer type */
          pmmdMol = (PMMD) pfbThis;
          c = (char) iIndex; 
          if (c == pmmdMol->pcMolName[0])
            {
               *ppdnmmFound = pmmdMol->pdnmmLink;
            }
      }
}


PDNMM LIBCALL MolFromPDBChain(PDNMS pdnmsThis,  char c)
{

    Int2 iTest;
    PDNMM pdnmmThis = NULL;
    
    if (!pdnmsThis)
      return NULL;
    
    iTest = TraverseMolecules( (DValNodePtr) pdnmsThis, 1, c,  
                          &pdnmmThis, FindMolByChar);
    return pdnmmThis;
}

   
CharPtr LIBCALL SeqStringFromMol(PDNMM pdnmmThis)
{
   PMMD pmmdThis = NULL; 
   PDNMG pdnmgThis = NULL;
   PMGD pmgdThis = NULL;
   Int4 iTest =0;
   CharPtr pcSeq = NULL;
   Int4 i = 0;
   
   if (!pdnmmThis)
     return NULL;
   
   pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
    
   if (pmmdThis->iResCount > 1) 
     {
       pcSeq = (CharPtr) MemNew((size_t) ((pmmdThis->iResCount + 1) * sizeof(Char)) );
       pdnmgThis = pmmdThis->pdnmgHead;
       while (pdnmgThis)
         {
	     pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
	     pcSeq[i] = pmgdThis->pcIUPAC[0];
	     i++;
	     pdnmgThis = pdnmgThis->next;
	 }
       pcSeq[i] = '\0';
     }
   return pcSeq;
}
 