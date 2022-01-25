#!/bin/csh -f
#
# Script to untar and make the NCBI toolkit on Solaris.
#
sed -f $NCBI/mk.csh.sed < $NCBI/ncbi.mk > ncbi.source.me
source ncbi.source.me
if ( "X" == "X$1") then
  echo USAGE $0 path to ncbi.tar file, including filename
  exit 1
endif

set tar_file = $1
set cwd = `pwd`

if ($tar_file != "-") then
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
   VIBFLAG="$NCBI_VIBFLAG"

make  SHELL="$NCBI_MAKE_SHELL" LCL="$NCBI_DEFAULT_LCL" RAN="$NCBI_RANLIB" CC="$NCBI_CC" \
   VIBLIBS="$NCBI_DISTVIBLIBS" \
   LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a \
   VIBFLAG="$NCBI_VIBFLAG"

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

if ($net_stat != 0 || $make_stat != 0 || $demo_stat != 0) then

   echo FAILURE primary make status = $make_stat, demo = $demo_stat, net = $net_stat
   exit 1
endif

