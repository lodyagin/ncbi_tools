# Microsoft Developer Studio Project File - Name="netentr" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=netentr - Win32 DebugMT
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "netentr.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "netentr.mak" CFG="netentr - Win32 DebugMT"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "netentr - Win32 DebugMT" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "netentr - Win32 DebugMT"

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
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /MTd /W3 /Gm /GR /Z7 /Od /I "..\..\..\..\.." /I "..\..\..\..\..\network\entrez\client" /I "..\..\..\..\..\object" /I "..\..\..\..\..\api" /I "..\..\..\..\..\asnlib" /I "..\..\..\..\..\corelib" /I "..\..\..\..\..\cdromlib" /I "..\..\..\..\..\asnstat" /I "..\..\..\..\..\biostruc" /I "..\..\..\..\..\network\nsclilib" /D "WIN32" /D "_DEBUG" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "netentr - Win32 DebugMT"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\..\network\entrez\client\netentr.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\entrez\client\netlib.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\entrez\client\objneten.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\..\network\entrez\client\netentr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\entrez\client\netlib.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\entrez\client\netpriv.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\entrez\client\objneten.h
# End Source File
# End Group
# End Target
# End Project
