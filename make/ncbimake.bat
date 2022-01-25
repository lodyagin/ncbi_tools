@echo off
rem  $Revision: 6.5 $
rem  +++++ Denis Vakatov, NCBI // vakatov@peony.nlm.nih.gov +++++

echo  *** Running:  %0  ***

if "%__PATH_16%"==""  set __PATH_16=c:\msvc
rem if "%__PATH_32%"==""  set __PATH_32=c:\msdev
if "%__PATH_32%"==""  set __PATH_32=D:\Program Files\DevStudio\SharedIDE\bin;D:\Program Files\DevStudio\VC
if "%__PATH_BOR%"==""  set __PATH_BOR=c:\bc45


echo  ***** Copy/Build/User flag  *****

if "%4"=="copy"  goto l_COPY
if "%4"==""      goto l_BUILD
if "%4"=="noexe" goto l_BUILD
if "%4"=="user"  goto l_USER
goto l_USAGE
:l_COPY
:l_BUILD
set __MAKEPATH=.
set __NCBIHOME=NNNNNNN
goto end_C_T
:l_USER
if "%6"==""  goto l_USAGE
if "%5"==""  goto l_USAGE
set __MAKEPATH=%6\make
set __NCBIHOME=NCBIHOME
:end_C_T


echo  ***** Debug/NoDebug  *****

if "%3"=="D"  goto l_DEB
if "%3"=="O"  goto l_OPT
goto l_USAGE
:l_DEB
set __DBUG=DBUG
goto end_D_O
:l_OPT
set __DBUG=DDDD
:end_D_O


echo  ***** Window/Console  *****

if "%2"=="W"  goto l_WIN
if "%2"=="C"  goto l_CON
goto l_USAGE
:l_WIN
set __CONSOLE=CCCCCCC
goto end_W_C
:l_CON
set __CONSOLE=CONSOLE
:end_W_C


echo  ***** MSVC-16/32/DLL32; Borland-16/32  *****

if "%1"=="16"   goto l_16
if "%1"=="32"   goto l_32
if "%1"=="32D"  goto l_32D
if "%1"=="BOR"  goto l_BOR
if "%1"=="B16"  goto l_B16

goto l_USAGE
:l_16
set __COMP=WWWWW
set BASE_PATH=%__PATH_16%
set __MAKE=nmake
goto end_COMP
:l_32
set __COMP=WIN32
set BASE_PATH=%__PATH_32%
set __MAKE=nmake
goto end_COMP
:l_32D
set __COMP=WIN32D
set BASE_PATH=%__PATH_32%
set __MAKE=nmake
goto end_COMP
:l_BOR
set __COMP=BOR
set BASE_PATH=%__PATH_BOR%
set __MAKE=make -N
goto end_COMP
:l_B16
set __COMP=B16
set BASE_PATH=%__PATH_BOR%
set __MAKE=make -N
:end_COMP


echo  ***** Set Environment  *****

set __PATH=%PATH%
set __LIB=%LIB%
set __BIN=%BIN%
set __INCLUDE=%INCLUDE%

set LIB=%BASE_PATH%\lib
set BIN=%BASE_PATH%\bin
set INCLUDE=%BASE_PATH%\include
set PATH=%BIN%;%PATH%


echo  ***** Run Makefile  *****

echo %__MAKE% -f %__MAKEPATH%\makefile.dos %__COMP%=1 %__CONSOLE%=1 %__DBUG%=1 MAKEUSER=%5.mak %__NCBIHOME%=%6 %4
%__MAKE% -f %__MAKEPATH%\makefile.dos %__COMP%=1 %__CONSOLE%=1 %__DBUG%=1 MAKEUSER=%5.mak %__NCBIHOME%=%6 %4


echo  *****  Reset Environment  *****

set LIB=%__LIB%
set BIN=%__BIN%
set INCLUDE=%__INCLUDE%
set PATH=%__PATH%

set __LIB=
set __BIN=
set __INCLUDE=
set __PATH=

set __DBUG=
set __COMP=
set __CONSOLE=
set __MAKE=
set __NCBIHOME=
set __MAKEPATH=

goto l_EXIT


echo  ***** Usage  *****

:l_USAGE
echo "Usage:"
echo " ncbimake {16|32|32D|BOR|B16} {W|C} {O|D} [noexe | copy | user <makefile> <toolkit_path>]"
echo " (16)/(32)-bit MSVC++;  (32D)LL MSVC++;  (BOR)land 32-bit;  (B16)Borland 16-bit"
echo " (W)indows/(C)onsole application"
echo " (D)ebug/(O)ptimized version"
echo " "
echo " Argument #4:"
echo "  (copy)  -- deploy all sources(must be run from MAKE directory)"
echo "  ()      -- build in the current directory(created by 'copy')"
echo "  (noexe) -- same as () but build libraries only(no executables)"
echo "  (user)  -- must be followed by(both):"
echo "     <makefile>      -- user's makefile basename (will be completed by .MAK)"
echo "     <toolkit_path>  -- path to the NCBI toolkit"
echo " "
echo "For more details please see in MAKE/README.DOS"


echo  ***** Exit  *****

:l_EXIT
set __PATH_16=
set __PATH_32=
set __PATH_BOR=

echo  *** Exiting:  %0  ***


