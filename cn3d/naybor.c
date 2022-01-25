/*   naybor.c
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
* File Name:  naybor.c
*
* Author:  Christopher Hogue
*
* Version Creation Date:   4/17/96
*
* $Revision: 6.2 $
*
* File Description: Cn3d file opening routines 
*                   
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* $Log: naybor.c,v $
* Revision 6.2  1998/08/26 18:28:39  kans
* fixed -v -fd warnings
*
* Revision 6.1  1998/03/06 23:19:25  lewisg
* codewarrior fixes
*
* Revision 6.0  1997/08/25 18:13:57  madden
* Revision changed to 6.0
*
* Revision 5.0  1996/05/28 14:05:44  ostell
* Set to revision 5.0
*
 * Revision 1.1  1996/04/18  18:02:59  epstein
 * Initial revision
 *
*
* ==========================================================================
*/

#include <vibrant.h>
#include <mmdbapi.h>
#include <naybor.h>
#include <cn3dmain.h>
#include <cn3dslct.h>
#include <cn3dwipe.h>
#include <cn3dxprt.h>
#include <cn3dsave.h>


static Boolean  Cn3D_Naybor_InUse = FALSE;
static IteM	Cn3D_iNeighbor;
static IteM     Cn3D_iNOpen;
static MenU	Cn3D_sNSave;
static MenU	Cn3D_sNExport;
static IteM	Cn3D_iNSelStruc;
static IteM	Cn3D_iNClearStruc;
static WindoW   Cn3D_wNOpen;

void LIBCALL Cn3D_EnableNayborOps(void)
{
 /*  Cn3D_DisableNayborOps();
   Enable(Cn3D_iNOpen);
   if (AreNeighborsOn()) 
   {
     if (GetSlaveModelstruc())
      {
        Enable(Cn3D_sNSave); 
        Enable(Cn3D_sNExport);
        Enable(Cn3D_iNSelStruc);
        Enable(Cn3D_iNClearStruc);
      }
   }
 */
}

void LIBCALL Cn3D_DisableNayborOps(void)
{      
/*	Disable(Cn3D_iNOpen);
  	Disable(Cn3D_sNSave); 
	Disable(Cn3D_sNExport); 
  	Disable(Cn3D_iNSelStruc);
	Disable(Cn3D_iNClearStruc);
 */	return;
}


static void Cn3D_NayborStatusProc(IteM i)
{
/*
  if (GetStatus(Cn3D_iNeighbor) == FALSE)
   {
     SetNeighborOff();
   }
  else
   {
     SetNeighborOn();
   }
  SetStatus(Cn3D_iNeighbor, AreNeighborsOn());  
  Cn3D_EnableNayborOps();
*/
}




static void Cn3D_OpenNayborDlg(IteM i)
{
    GrouP   g;
    CharPtr name = NULL;
    Char    buf[20];
    
    buf [0] = '\0';
    if (Cn3D_Naybor_InUse) return;
    else Cn3D_Naybor_InUse = TRUE;
    name = StringSave(GetStrucStrings(GetMasterModelstruc(), PDB_ACC));
    StringCpy(buf,"Neighbors to ");
    if (!name) name = StringSave("UNKN");
    StringNCat(buf, name, 4);
    MemFree(name);
    Cn3D_wNOpen = FixedWindow(-30, -20,  -10,  -10, "Open Neighbor",  NULL);
    g = NormalGroup(Cn3D_wNOpen, 2, 1, buf ,  systemFont, NULL);

    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 20);  
   
    Show(Cn3D_wNOpen);
    return;
}
/*
BiostrucAnnotSetPtr CDECL NetEntrezBiostrucAnnotSetGet (DocUid gi)



After the libraries are rebuilt (deferred until later today because Jim
is in the midst of a lot of changes), the CD-ROM-based Entrez API will support
structure alignments via a kludge in cdentrez.c which gets these alignments
from Hitomi at:
  /net/clobber/usr/people6/ohkawa/mdb/biostruc/VASTDATA/ *.align

Here's what's there at present:

gold> ls /net/clobber/usr/people6/ohkawa/mdb/biostruc/VASTDATA/ *.align
/net/clobber/usr/people6/ohkawa/mdb/biostruc/VASTDATA/1AAZ.align
/net/clobber/usr/people6/ohkawa/mdb/biostruc/VASTDATA/1ABA.align
/net/clobber/usr/people6/ohkawa/mdb/biostruc/VASTDATA/1MCP.align
/net/clobber/usr/people6/ohkawa/mdb/biostruc/VASTDATA/1PPI.align
/net/clobber/usr/people6/ohkawa/mdb/biostruc/VASTDATA/2MB5.align

Greg:
When you get back, please modify cdentrez.c so that it gets the
alignments from the 'right place' (the disk-based "CD-ROM" dataset).
I made the assumption for now that you'll be using the "extra" field
of the structure Docsum record to store an offset to the alignment
record ... see "rcsdiff cdentrez.c" for details.

The Network Entrez API also supports this, but I want to defer installing
a new Network Entrez server until after Greg has replaced my kludge with
the real implementation.


I just installed a new Network Entrez server for which the alignments are
available via the EntrezBiostrucAnnotSetGet() function call.  You can test
this with MMDB ID 164.  New alignments will be available in this fashion
after they are delivered to Greg.

162
164
1649
2788
3575

*/

MenU LIBCALL Cn3D_NayborSub (MenU m)
{
  MenU s;
  s = SubMenu (m, "Neighbor");
  Cn3D_iNeighbor  = StatusItem (s, " Mode/M", NULL /*Cn3D_NayborStatusProc */);
  SeparatorItem(s);
  Cn3D_iNOpen = CommandItem (s, "Open/O", NULL /* Cn3D_OpenNayborDlg */);
  SeparatorItem(s);
  Cn3D_iNSelStruc = CommandItem(s, "Active /A",NULL /* Cn3D_SelectDlg */);
       	/* see cn3dslct.c */
  Cn3D_iNClearStruc = CommandItem(s, "Clear /C", NULL /* Cn3D_ClearSelProc */);
	/* see cn3dwipe	*/

  SeparatorItem(s);
  Cn3D_sNSave = Cn3D_SaveSub (s);
  Cn3D_sNExport = Cn3D_ExportSub (s);
  Cn3D_EnableNayborOps();

  return s;

}
