#
# $Id: r6k.ncbi.mk,v 1.1 2000/01/11 17:03:10 beloslyu Exp $
#
# As per ddas@us.ibm.com 08/19/1999
#
NCBI_DEFAULT_LCL = r6k
NCBI_MAKE_SHELL = /bin/sh
NCBI_CC = xlc
NCBI_CFLAGS1 = -c -Xa
NCBI_LDFLAGS1 = -Xa
NCBI_OPTFLAG = -O
NCBI_INCDIR = /usr/ncbi/include
NCBI_LIBDIR = /usr/ncbi/lib
NCBI_ALTLIB = /usr/ncbi/altlib
#will work only when you have Motif installed!
NCBI_VIBFLAG = -I/usr/X11R6/include -L/usr/X11R6/lib -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXmu -lXt -lX11
NCBI_DISTVIBLIBS = -L/usr/X11R6/lib -lXm -lXmu -lXt -lX11
NCBI_OTHERLIBS = -lm
NCBI_RANLIB = ranlib
# Used by makedis.csh
NCBI_MT_OTHERLIBS = -lpthread
NCBI_THREAD_OBJ = ncbithr.o
NETENTREZVERSION = 2.02c2ASN1SPEC6 

