# Microsoft Developer Studio Project File - Name="ncbicn3d" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ncbicn3d - Win32 DebugMT
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ncbicn3d.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ncbicn3d.mak" CFG="ncbicn3d - Win32 DebugMT"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ncbicn3d - Win32 DebugMT" (based on "Win32 (x86) Static Library")
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
# ADD CPP /nologo /MTd /W3 /GR /Zi /Od /I "..\..\..\.." /I "..\..\..\..\corelib" /I "..\..\..\..\cn3d" /I "..\..\..\..\vibrant" /I "..\..\..\..\object" /I "..\..\..\..\biostruc" /I "..\..\..\..\asnlib" /I "..\..\..\..\asnstat" /I "..\..\..\..\desktop" /I "..\..\..\..\access" /I "..\..\..\..\api" /I "..\..\..\..\ddv" /I "..\..\..\..\tools" /I "..\..\..\..\cdromlib" /D "WIN32" /D "_DEBUG" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
# Begin Target

# Name "ncbicn3d - Win32 DebugMT"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\cn3d\algorend.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dentr.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmatn.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmesh.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmodl.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmsel.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmsg.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dopen.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dsave.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dshim.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dstyl.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dwin.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dxprt.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\seqcons.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\cn3d\algorend.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dentr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmain.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmesh.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmodl.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmsel.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dmsg.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dopen.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dsave.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dshim.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dstyl.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\cn3dxprt.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\resource.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\cn3d\seqcons.h
# End Source File
# End Group
# End Target
# End Project
