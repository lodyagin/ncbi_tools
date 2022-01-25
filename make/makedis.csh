#!/bin/csh -f
#
# $Id: makedis.csh,v 1.20 1999/03/18 17:33:49 beloslyu Exp $
#
##                            PUBLIC DOMAIN NOTICE                          
#               National Center for Biotechnology Information
#                                                                          
#  This software/database is a "United States Government Work" under the   
#  terms of the United States Copyright Act.  It was written as part of    
#  the author's official duties as a United States Government employee and 
#  thus cannot be copyrighted.  This software/database is freely available 
#  to the public for use. The National Library of Medicine and the U.S.    
#  Government have not placed any restriction on its use or reproduction.  
#                                                                          
#  Although all reasonable efforts have been taken to ensure the accuracy  
#  and reliability of the software and data, the NLM and the U.S.          
#  Government do not and cannot warrant the performance or results that    
#  may be obtained by using this software or data. The NLM and the U.S.    
#  Government disclaim all warranties, express or implied, including       
#  warranties of performance, merchantability or fitness for any particular
#  purpose.                                                                
#                                                                          
#  Please cite the author in any work or product based on this material.   
# Author: Karl Sirotkin <sirotkin@ncbi.nlm.nih.gov>
#
#
# Script to untar and make the NCBI toolkit on Solaris.
#

set MFLG=""

if ("X$1" == "X-n") then
	set MFLG="-n"
	shift	
endif

set tar_file = $1
set cwd = `pwd`

# do we need to extract the tar archive?
if ( "X$tar_file" != "X" && "$tar_file" != "-") then
	if (! -r "$tar_file") then
		echo Unable to find the file "$tar_file"
		exit 1
	endif

	if (-d "ncbi") then
		echo "ncbi directory already exists, please remove or rename it"
		exit 2
	endif

	ls -l $tar_file
	tar xvf $tar_file
else
	# make sure that ncbi/build directory exists
	if ( ! -d "ncbi/build" ) then
		echo 'ncbi/build directory should exist. Did you extract ncbi.tar.Z?'
		exit 2
	endif
endif

set os=`uname -s`
switch ($os)
case SunOS:
	switch (`uname -r`)
	case "4.1*":
		set platform=sun
		breaksw
	default:
		if ( `uname -p` == i386 ) then
			set platform=solarisintel
		else
			set platform=solaris
		endif
		breaksw
	endsw
	breaksw
case IRIX*:
	switch (`uname -r`)
	case "5.*":
		set platform=sgi5
		breaksw
	case "4.*":
		set platform=sgi4
		breaksw
	case "6.*":
		set platform=sgi
		breaksw
	default:
		set platform=sgi
		breaksw
	endsw
	breaksw
case OSF1:
	set platform=alphaOSF1
	breaksw
case Linux:
	set platform=linux
	breaksw
default:
	echo Platform not found : `uname -a`
	goto BADPLATFORM
	breaksw
endsw

echo platform is $platform

set NCBI_DOT_MK = ncbi/platform/${platform}.ncbi.mk

if (! -r "$NCBI_DOT_MK") then
  goto BADPLATFORM
endif

set noglob
# take the file $NCBI_DOT_MK and convert it to be suitable for csh eval:
# (1) remove comments at the beginning of the lines
# (2) change variable referenses to be in curly brackets - $AAA -> ${AAA}
# (3) remove excessive spaces around the equal sign
# (4) change Makefile assignments to csh ones: AAA=bb cc -> set AAA = "bb cc"
eval `sed -e 's/^ *#.*//g' -e 's/\$(\([a-zA-Z_]*\))/\${\1}/g' -e 's/ *= */=/g' -e 's/^\([^=]*\)=\(.*\)$/set \1 = "\2";/' < $NCBI_DOT_MK`
unset noglob

cd ncbi/build
cp ../make/*.unx .
mv makeall.unx makefile

#  Inherited to this system is the requirement to use:
#    TO USE VIBRANT
# to have for makeall, this line
#  LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a 
#    and for the makenet, this symbol
#  BLIB31=libvibnet.a 
#
set make="make"
if ("$platform" == "solaris" || "$platform" == "solarisintel") then
	set tmp=`/usr/bin/which dmake|sed -e 's/^\(...\).*/\1/'`
	if ("$tmp" != "no ") then
		set make="dmake -j 2"
	endif
endif

# it's not working reliably
#if ("$platform" == "alphaOSF1") then
#	set make="make -j 2"
#endif

set CMD='$make $MFLG SHELL=\"$NCBI_MAKE_SHELL\" LCL=\"$NCBI_DEFAULT_LCL\" \
   RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\" VIBLIBS=\"$NCBI_DISTVIBLIBS\" \
   LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a \
   LIB20=libncbidesk.a VIBFLAG=\"$NCBI_VIBFLAG\" all'
eval echo $CMD
eval echo $CMD | sh 

set make_stat = $status

set CMD='$make $MFLG -f makedemo.unx SHELL=\"$NCBI_MAKE_SHELL\" \
   LCL=\"$NCBI_DEFAULT_LCL\" RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\" \
   VIBLIBS=\"$NCBI_DISTVIBLIBS\" LIB4=-lvibrant VIBFLAG=\"$NCBI_VIBFLAG\"'
eval echo $CMD
eval echo $CMD | sh 

set demo_stat = $status

rm -f blastall blastpgp seedtop

#
# In case platform supports multi-threading, remake the apps which
# should be multithreaded, if at all possible.
#  Might repeat what is done above on some platforms.
#

set CMD='$make $MFLG -f makedemo.unx SHELL=\"$NCBI_MAKE_SHELL\" \
   LCL=\"$NCBI_DEFAULT_LCL\" RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\"  \
   VIBLIBS=\"$NCBI_DISTVIBLIBS\" LIB4=-lvibrant \
   THREAD_OBJ=$NCBI_THREAD_OBJ THREAD_OTHERLIBS=$NCBI_MT_OTHERLIBS \
   VIBFLAG=\"$NCBI_VIBFLAG\" blastall blastpgp seedtop'
eval echo $CMD
eval echo $CMD | sh 

set threaded_demo_stat = $status

set CMD='$make $MFLG -f makenet.unx SHELL=\"$NCBI_MAKE_SHELL\" \
   CC=\"$NCBI_CC\" RAN=\"$NCBI_RANLIB\" VIBLIBS=\"$NCBI_DISTVIBLIBS\" \
   VIBFLAG=\"$NCBI_VIBFLAG\" \
   VIB=\"Psequin Nentrez Cn3D powblast pblcmd blastcl3\" \
   BLIB31=libvibnet.a OTHERLIBS=\"$NCBI_OTHERLIBS\"'
eval echo $CMD
eval echo $CMD | sh 

set net_stat = $status

if ($net_stat != 0 || $make_stat != 0 || $demo_stat != 0 || $threaded_demo_stat != 0) then

   echo FAILURE primary make status = $make_stat, demo = $demo_stat, net = $net_stat, threaded_demo_stat = $threaded_demo_stat
   exit 1
else
   # we are in ncbi/build directory now. Let us make the VERSION file
   echo putting date stamp to the file ../VERSION
   date > ../VERSION
   exit 0
endif

BADPLATFORM:
  echo 'Your platform is not supported.'
  echo 'To port ncbi toolkit to your platform consult'
  echo 'the files platform/*.ncbi.mk'

exit 0
