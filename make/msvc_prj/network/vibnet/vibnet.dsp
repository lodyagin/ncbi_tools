# Microsoft Developer Studio Project File - Name="vibnet" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=vibnet - Win32 DebugMT
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "vibnet.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "vibnet.mak" CFG="vibnet - Win32 DebugMT"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "vibnet - Win32 DebugMT" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe
# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "DebugMT"
# PROP BASE Intermediate_Dir "DebugMT"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "DebugMT"
# PROP Intermediate_Dir "DebugMT"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /MTd /W3 /GR /Zi /Od /I "..\..\..\.." /I "..\..\..\..\vibrant" /I "..\..\..\..\asnlib" /I "..\..\..\..\api" /I "..\..\..\..\asnstat" /I "..\..\..\..\corelib" /I "..\..\..\..\desktop" /I "..\..\..\..\connect" /I "..\..\..\..\cdromlib" /I "..\..\..\..\object" /I "..\..\..\..\biostruc" /I "..\..\..\..\network\vibnet" /I "..\..\..\..\cn3d" /I "..\..\..\..\tools" /I "..\..\..\..\ddv" /D "WIN32" /D "_DEBUG" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
# Begin Target

# Name "vibnet - Win32 DebugMT"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\network\vibnet\docsum.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\network\vibnet\netcnfg.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\network\vibnet\trmlst.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\network\vibnet\entrez.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\network\vibnet\netcnfg.h
# End Source File
# End Group
# End Target
# End Project
