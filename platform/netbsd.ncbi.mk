#
# $Id: netbsd.ncbi.mk,v 1.2 1999/08/12 18:32:41 beloslyu Exp $
#
# That platform is not supported and is contributed to NCBI by 
#                     Brook Milligan <brook@trillium.NMSU.Edu>
#
NCBI_DEFAULT_LCL = bsd
NCBI_MAKE_SHELL = /bin/sh
NCBI_CC = gcc -pipe
NCBI_CFLAGS1 = -c
NCBI_LDFLAGS1 = -O
NCBI_OPTFLAG = -O
NCBI_BIN_MASTER = /home/coremake/ncbi/bin
NCBI_BIN_COPY = /home/coremake/ncbi/bin
NCBI_INCDIR = /home/coremake/ncbi/include
NCBI_LIBDIR = /home/coremake/ncbi/lib
NCBI_ALTLIB = /home/coremake/ncbi/altlib
#will work only when you have Motif installed!
NCBI_VIBFLAG = -I/usr/X11R6/include -L/usr/X11R6/lib -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXmu -lXp -lXpm -lXt -lX11 -lXext 
NCBI_DISTVIBLIBS = -L/usr/X11R6/lib /usr/X11R6/lib/libXm.a /usr/X11R6/lib/libXp.a /usr/X11R6/lib/libXpm.a -lXmu -lXt -lX11 -lXext
NCBI_OTHERLIBS = -lm
NCBI_RANLIB = ranlib
# Used by makedis.csh
NCBI_MT_OTHERLIBS = 
NCBI_THREAD_OBJ = ncbithr.o
