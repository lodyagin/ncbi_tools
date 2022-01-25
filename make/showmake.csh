#!/bin/csh -f
#
# $Id: showmake.csh,v 1.4 1999/06/07 18:13:48 beloslyu Exp $
#
#
# Show how to make the NCBI toolkit 
#
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
case NetBSD:
	set platform=netbsd
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

echo 'cd ncbi/build'
echo 'cp ../make/*.unx .'
echo 'mv makeall.unx makefile'

echo make  SHELL=\"$NCBI_MAKE_SHELL\" LCL=\"$NCBI_DEFAULT_LCL\" RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\" \
   VIBLIBS=\"$NCBI_DISTVIBLIBS\" \
   LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a \
   VIBFLAG=\"$NCBI_VIBFLAG\" all
echo ''
echo make -f makedemo.unx SHELL=\"$NCBI_MAKE_SHELL\" LCL=\"$NCBI_DEFAULT_LCL\" \
   RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\"  \
   VIBLIBS=\"$NCBI_DISTVIBLIBS\" \
   LIB4=-lvibrant \
   VIBFLAG=\"$NCBI_VIBFLAG\"
echo ''
echo make -f makedemo.unx SHELL=\"$NCBI_MAKE_SHELL\" LCL=\"$NCBI_DEFAULT_LCL\" \
   RAN=\"$NCBI_RANLIB\" CC=\"$NCBI_CC\"  \
   VIBLIBS=\"$NCBI_DISTVIBLIBS\" \
   LIB4=-lvibrant \
   THREAD_OBJ=$NCBI_THREAD_OBJ THREAD_OTHERLIBS=$NCBI_MT_OTHERLIBS \
   VIBFLAG=\"$NCBI_VIBFLAG\" blastall blastpgp seedtop
echo ''
echo make -f makenet.unx SHELL=\"$NCBI_MAKE_SHELL\" CC=\"$NCBI_CC\" \
   RAN=\"$NCBI_RANLIB\" \
   VIBLIBS=\"$NCBI_DISTVIBLIBS\" \
   VIBFLAG=\"$NCBI_VIBFLAG\" \
   VIB=\"Psequin Nentrez Cn3D blastcl3\" BLIB31=libvibnet.a OTHERLIBS=\"$NCBI_OTHERLIBS\"

exit 0

BADPLATFORM:
echo 'Your platform is not supported.'
echo 'To port ncbi toolkit to your platform consult'
echo 'the files platform/*.ncbi.mk'

exit 1
