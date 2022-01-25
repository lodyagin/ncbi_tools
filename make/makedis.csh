#!/bin/csh -f
#
# $Id: makedis.csh,v 1.13 1999/01/25 18:29:02 beloslyu Exp $
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
# Author: Karl Sirotkin

#
# Script to untar and make the NCBI toolkit on Solaris.
#
if (! $?NCBI) then
  goto BEGINNER
endif
if (! -r $NCBI/ncbi.mk) then
  goto BEGINNER
endif
sed -f $NCBI/mk.csh.sed < $NCBI/ncbi.mk > ncbi.source.me
source ncbi.source.me
if ( "X" == "X$1") then
  echo USAGE $0 path to ncbi.tar file, including filename
  exit 1
endif

set tar_file = $1
set cwd = `pwd`

if ( X$tar_file != "X-") then
	if (! -e $tar_file) then
		echo Unable to find $tar_file
		exit 1
	endif

	if (-d "ncbi") then
		echo "ncbi directory already exists, please remove or rename"
		exit 2
	endif

 \ls -l $tar_file
	tar -xvf $tar_file

endif

cd ncbi/build
cp ../make/*.unx .
mv makeall.unx makefile
#
#  Inherited to this system is the requirement to use:
#    TO USE VIBRANT
# to have for makeall, this line
#  LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a 
#    and for the makenet, this symbol
#  BLIB31=libvibnet.a 
#
echo make  SHELL="$NCBI_MAKE_SHELL" LCL="$NCBI_DEFAULT_LCL" RAN="$NCBI_RANLIB" CC="$NCBI_CC" \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a \
   VIBFLAG="$NCBI_VIBFLAG" all

make  SHELL="$NCBI_MAKE_SHELL" LCL="$NCBI_DEFAULT_LCL" RAN="$NCBI_RANLIB" CC="$NCBI_CC" \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a \
   VIBFLAG="$NCBI_VIBFLAG" all

set make_stat = $status

echo make -f makedemo.unx SHELL="$NCBI_MAKE_SHELL" LCL="$NCBI_DEFAULT_LCL" \
   RAN="$NCBI_RANLIB" CC="$NCBI_CC"  \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   LIB4=-lvibrant \
   VIBFLAG="$NCBI_VIBFLAG"

make -f makedemo.unx SHELL="$NCBI_MAKE_SHELL" LCL="$NCBI_DEFAULT_LCL" \
   RAN="$NCBI_RANLIB" CC="$NCBI_CC"  \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   LIB4=-lvibrant \
   VIBFLAG="$NCBI_VIBFLAG"

set demo_stat = $status

rm -f blastall blastpgp seedtop

#
# In case platform supports multi-threading, remake the apps which
# should be multithreaded, if at all possible.
#  Might repeat what is done above on some platforms.
#

echo make -f makedemo.unx SHELL="$NCBI_MAKE_SHELL" LCL="$NCBI_DEFAULT_LCL" \
   RAN="$NCBI_RANLIB" CC="$NCBI_CC"  \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   LIB4=-lvibrant \
   THREAD_OBJ=$NCBI_THREAD_OBJ THREAD_OTHERLIBS=$NCBI_MT_OTHERLIBS \
   VIBFLAG="$NCBI_VIBFLAG" blastall blastpgp seedtop

make -f makedemo.unx SHELL="$NCBI_MAKE_SHELL" LCL="$NCBI_DEFAULT_LCL" \
   RAN="$NCBI_RANLIB" CC="$NCBI_CC"  \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   LIB4=-lvibrant \
   THREAD_OBJ=$NCBI_THREAD_OBJ THREAD_OTHERLIBS=$NCBI_MT_OTHERLIBS \
   VIBFLAG="$NCBI_VIBFLAG" blastall blastpgp seedtop

set threaded_demo_stat = $status

echo make -f makenet.unx SHELL="$NCBI_MAKE_SHELL" CC="$NCBI_CC" \
   RAN="$NCBI_RANLIB" \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   VIBFLAG="$NCBI_VIBFLAG" \
   VIB="Psequin Nentrez Cn3D blastcl3" BLIB31=libvibnet.a OTHERLIBS="$NCBI_OTHERLIBS"

make -f makenet.unx SHELL="$NCBI_MAKE_SHELL" CC="$NCBI_CC" \
   RAN="$NCBI_RANLIB" \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   VIBFLAG="$NCBI_VIBFLAG" \
   VIB="Psequin Nentrez Cn3D blastcl3" BLIB31=libvibnet.a OTHERLIBS="$NCBI_OTHERLIBS"

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

BEGINNER:
  echo please read makedis.csh which will either be in the make
  echo subdirectory, or available from ftp from
  echo ncbi.nlm.nih.gov at toolbox/readme.unx
  echo it has information on how to optain the ncbi.mk file
  echo other helpful hints.

exit 0
