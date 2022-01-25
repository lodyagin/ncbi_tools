#!/bin/csh -f
#
# $Id: showmake.csh,v 1.2 1999/01/08 18:04:57 coremake Exp $
#
#
# Show how to make the NCBI toolkit 
#
sed -f $NCBI/mk.csh.sed < $NCBI/ncbi.mk > ncbi.source.me
source ncbi.source.me
echo 'cd ncbi/build'
echo 'cp ../make/*.unx .'
echo mv makeall.unx makefile
echo make  SHELL=\"$NCBI_MAKE_SHELL\" LCL=\"$NCBI_DEFAULT_LCL\" RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\"  VIBLIBS=\"$NCBI_DISTVIBLIBS\"  LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a  VIBFLAG=\"$NCBI_VIBFLAG\"


echo make -f makedemo.unx SHELL=\"$NCBI_MAKE_SHELL\" LCL=\"$NCBI_DEFAULT_LCL\"  RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\"   VIBLIBS=\"$NCBI_DISTVIBLIBS\"  LIB4=-lvibrant  VIBFLAG=\"$NCBI_VIBFLAG\"

echo make -f makenet.unx SHELL=\"$NCBI_MAKE_SHELL\" CC=\"$NCBI_CC\"  VIBLIBS=\"$NCBI_DISTVIBLIBS\"  VIBFLAG=\"$NCBI_VIBFLAG\"  VIB=\"Psequin Nentrez Cn3D\" BLIB31=\"libvibnet.a\" OTHERLIBS=\"$NCBI_OTHERLIBS\"
