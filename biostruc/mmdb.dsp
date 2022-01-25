# Microsoft Developer Studio Project File - Name="mmdb" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=mmdb - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mmdb.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mmdb.mak" CFG="mmdb - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mmdb - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "mmdb - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "mmdb - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "mmdb - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I ".." /I "..\ddv" /I "..\cn3d" /I "..\access" /I "..\asnstat" /I "..\connect\lbapi" /I "..\connect" /I "..\asnlib" /I "..\vibrant" /I "..\biostruc" /I "..\object" /I "..\api" /I "..\cdromlib" /I "..\desktop" /I "..\tools" /I "..\corelib" /I "..\network\taxon1\common" /I "..\network\vibnet" /I "..\network\entrez\client" /I "..\network\nsclilib" /I "..\network\medarch\client" /I "..\network\id1arch" /I "..\network\taxon1\taxon2" /I "..\network\blast3\client" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Debug\ncbimmdb.lib"

!ENDIF 

# Begin Target

# Name "mmdb - Win32 Release"
# Name "mmdb - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\corematx.c
# End Source File
# Begin Source File

SOURCE=.\dvncode.c
# End Source File
# Begin Source File

SOURCE=.\mmdbapi.c
# End Source File
# Begin Source File

SOURCE=.\mmdbapi1.c
# End Source File
# Begin Source File

SOURCE=.\mmdbapi2.c
# End Source File
# Begin Source File

SOURCE=.\mmdbapi3.c
# End Source File
# Begin Source File

SOURCE=.\mmdbapi4.c
# End Source File
# Begin Source File

SOURCE=.\mmdbentr.c
# End Source File
# Begin Source File

SOURCE=.\objmmdb1.c
# End Source File
# Begin Source File

SOURCE=.\objmmdb2.c
# End Source File
# Begin Source File

SOURCE=.\objmmdb3.c
# End Source File
# Begin Source File

SOURCE=.\prunebsc.c
# End Source File
# Begin Source File

SOURCE=.\vastsubs.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\corematx.h
# End Source File
# Begin Source File

SOURCE=.\dvncode.h
# End Source File
# Begin Source File

SOURCE=.\mmdbapi.h
# End Source File
# Begin Source File

SOURCE=.\mmdbapi1.h
# End Source File
# Begin Source File

SOURCE=.\mmdbapi2.h
# End Source File
# Begin Source File

SOURCE=.\mmdbapi3.h
# End Source File
# Begin Source File

SOURCE=.\mmdbapi4.h
# End Source File
# Begin Source File

SOURCE=.\mmdbdata.h
# End Source File
# Begin Source File

SOURCE=.\objmmdb1.h
# End Source File
# Begin Source File

SOURCE=.\objmmdb2.h
# End Source File
# Begin Source File

SOURCE=.\objmmdb3.h
# End Source File
# Begin Source File

SOURCE=.\prunebsc.h
# End Source File
# Begin Source File

SOURCE=.\strimprt.h
# End Source File
# End Group
# End Target
# End Project
