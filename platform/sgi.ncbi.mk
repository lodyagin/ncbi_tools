#
# $Id: sgi.ncbi.mk,v 1.6 1999/08/12 18:32:41 beloslyu Exp $
#
NCBI_MAKE_SHELL = /bin/sh
NCBI_DEFAULT_LCL = sgi
NCBI_CC = cc -mips1 `uname -r | sed 's/6\..*/-32/;s/5\..*//'` 
NCBI_SYBASE = /usr/people/sybase_10.0.3
NCBI_SYBLIBS = -L$(NCBI_SYBASE)/lib -lsybdb
NCBI_SYBLIBS_STATIC = $(NCBI_SYBASE)/lib/libsybdb.a
NCBI_SYBFLAG = -I$(NCBI_SYBASE)/include -L$(NCBI_SYBASE)/lib
NCBI_CFLAGS1 = -c
NCBI_LDFLAGS1 = -O
NCBI_OPTFLAG = -O
NCBI_INCDIR = /usr/ncbi/ncbi/include/NCBI
NCBI_LIBDIR = /usr/ncbi/ncbi/lib
NCBI_ALTLIB =  /usr/ncbi/ncbi/altlib
NCBI_VIBFLAG = -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXmu -lXt -lX11
#
#if you want to link Motif statically, use that line
#NCBI_DISTVIBLIBS = -B static -lXm -B dynamic -lXmu -lXt -lX11 -lgen
#but if you want to link Motif dynamically, use that line
NCBI_DISTVIBLIBS = -lXm -lXmu -lXt -lX11 -lgen
#
NCBI_OTHERLIBS = -lm -lPW -lrpcsvc
NCBI_OTHERLIBS_MT = $(NCBI_OTHERLIBS) -lpthread
# NCBI_MT_OTHERLIBS & NCBI_THREAD_OBJ are only used by master makefiles
NCBI_MT_OTHERLIBS = -lpthread
NCBI_THREAD_OBJ = ncbithr.o
NCBI_THR_OBJ = $(NCBI_LIBDIR)/ncbithr.o
NCBI_THR_ALTOBJ = $(NCBI_ALTLIB)/ncbithr.o

NCBI_RANLIB = ls -l

NCBI_OPTIONAL_LIBS = BLIB42=libctutils.a BLIB44=libidload.a

NCBI_ALTSRC =  /usr/ncbi/ncbi/altsrc

SEQUINDOC = /am/ncbiapdata/sequin.htm
SEQUINIMAGES = /am/ncbiapdata/images

NETENTREZVERSION = 2.02c2ASN1SPEC6 
