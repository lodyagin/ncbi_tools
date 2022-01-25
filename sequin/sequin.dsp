# Microsoft Developer Studio Project File - Name="sequin" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=sequin - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "sequin.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "sequin.mak" CFG="sequin - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "sequin - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "sequin - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "sequin - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386

!ELSEIF  "$(CFG)" == "sequin - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I ".." /I "..\ddv" /I "..\cn3d" /I "..\access" /I "..\asnstat" /I "..\connect\lbapi" /I "..\connect" /I "..\asnlib" /I "..\vibrant" /I "..\biostruc" /I "..\object" /I "..\api" /I "..\cdromlib" /I "..\desktop" /I "..\tools" /I "..\corelib" /I "..\network\taxon1\common" /I "..\network\vibnet" /I "..\network\entrez\client" /I "..\network\nsclilib" /I "..\network\medarch\client" /I "..\network\id1arch" /I "..\network\taxon1\taxon2" /I "..\network\blast3\client" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_NETENT_" /D "_OPENGL" /FR /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib opengl32.lib glu32.lib glaux.lib wsock32.lib vibnet.lib ncbidesk.lib ncbimmdb.lib ncbitool.lib ncbinacc.lib netentr.lib netcli.lib vibrant.lib ncbicdr.lib ncbiobj.lib ncbi.lib ncbicn3d.lib ncbitool.lib ncbitxc2.lib ncbiid1.lib medarch.lib netblast.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept /libpath:"..\cdromlib\debug" /libpath:"..\cdromlib\ncbicdr___Win32_Debug" /libpath:"..\corelib\debug" /libpath:"..\vibrant\debug" /libpath:"..\object\debug" /libpath:"..\biostruc\debug" /libpath:"..\desktop\debug" /libpath:"..\tools\debug" /libpath:"..\network\entrez\client\debug" /libpath:"..\network\vibnet\debug" /libpath:"..\network\nsclilib\debug" /libpath:"..\network\taxon1\taxon2\debug" /libpath:"..\network\taxon1\common\debug" /libpath:"..\network\id1arch\debug" /libpath:"..\network\medarch\client\debug" /libpath:"..\network\blast3\client\debug" /libpath:"..\cn3d\ncbicn3d___Win32_Debug"
# SUBTRACT LINK32 /pdb:none

!ENDIF 

# Begin Target

# Name "sequin - Win32 Release"
# Name "sequin - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\link\mswin\sequin.rc
# End Source File
# Begin Source File

SOURCE=.\sequin1.c
# End Source File
# Begin Source File

SOURCE=.\sequin2.c
# End Source File
# Begin Source File

SOURCE=.\sequin3.c
# End Source File
# Begin Source File

SOURCE=.\sequin4.c
# End Source File
# Begin Source File

SOURCE=.\sequin5.c
# End Source File
# Begin Source File

SOURCE=.\sequin6.c
# End Source File
# Begin Source File

SOURCE=.\sequin7.c
# End Source File
# Begin Source File

SOURCE=.\sequin8.c
# End Source File
# Begin Source File

SOURCE=.\sequin9.c
# End Source File
# Begin Source File

SOURCE=.\sequinx.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\sequin.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
