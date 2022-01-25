# Microsoft Developer Studio Project File - Name="regexp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=regexp - Win32 DebugMT
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "regexp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "regexp.mak" CFG="regexp - Win32 DebugMT"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "regexp - Win32 DebugMT" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "regexp - Win32 DebugMT"

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
# ADD CPP /nologo /MTd /W3 /Gm /GR /Z7 /Od /I "..\..\.." /I "..\..\..\corelib" /D "WIN32" /D "_DEBUG" /YX /FD /GZ /c
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

# Name "regexp - Win32 DebugMT"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\regexp\chartables.c
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\get.c
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\maketables.c
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\pcre.c
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\pcreposix.c
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\study.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\regexp\pcre.h
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\pcre_config.h
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\pcre_internal.h
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\pcreposix.h
# End Source File
# End Group
# Begin Group "Doc"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\regexp\doc\AUTHORS
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\doc\LICENCE
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\doc\NEWS
# End Source File
# Begin Source File

SOURCE="..\..\..\regexp\doc\NON-UNIX-USE"
# End Source File
# Begin Source File

SOURCE=..\..\..\regexp\doc\README
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\regexp\README
# End Source File
# End Target
# End Project
