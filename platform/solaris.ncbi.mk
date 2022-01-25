#
# $Id: solaris.ncbi.mk,v 1.3 1999/04/21 13:21:43 beloslyu Exp $
#
NCBI_DEFAULT_LCL = sol
NCBI_MAKE_SHELL = /usr/bin/sh
NCBI_CC = cc -xildoff
NCBI_CFLAGS1 = -c -Xa
NCBI_LDFLAGS1 = -Xa
NCBI_OPTFLAG = -O -fast
# following 2 lines are temporary; J. Epstein, 8/11/97
NCBI_INCDIR = /netopt/ncbi_tools/ver0.0/ncbi/include
NCBI_LIBDIR = /netopt/ncbi_tools/ver0.0/ncbi/lib
NCBI_ALTLIB = /netopt/ncbi_tools/ver0.0/ncbi/altlib
NCBI_SHLIB = /netopt/ncbi_tools/ver0.0/ncbi/shlib -R/netopt/ncbi_tools/ver0.0/ncbi/shlib
NCBI_CLLIB = /netopt/ncbi_tools/ver0.0/ncbi/cllib
NCBI_CLINCDIR = ./ -I/netopt/ncbi_tools/ver0.0/ncbi/cllib -I/netopt/ncbi_tools/ver0.0/ncbi/include
NCBI_ALTSRC = /netopt/ncbi_tools/ver0.0/ncbi/altsrc
#PRE-CDE NCBI_VIBFLAG = -I/usr/openwin/include -I/netopt/SUNWmotif/include -L/usr/openwin/lib -L/netopt/SUNWmotif/lib -DWIN_MOTIF
NCBI_VIBFLAG = -I/usr/openwin/include -I/usr/dt/include -L/usr/openwin/lib -L/usr/dt/lib -DWIN_MOTIF
NCBI_VIBLIBS = -R/usr/dt/lib -R/usr/openwin/lib -lXm -lXmu -lXt -lX11 -lXext
#PRE-CDE NCBI_DISTVIBLIBS = -Bstatic -lXm -Bdynamic -lXmu -lXt -lX11 -lXext
NCBI_DISTVIBLIBS = -R/usr/dt/lib:/usr/openwin/lib -lXm -lXmu -lXt -lX11 -lXext
NCBI_OTHERLIBS = -lresolv -lsocket -lrpcsvc -lnsl -lgen -lm
# NCBI_MT_OTHERLIBS & NCBI_THREAD_OBJ are only used by master makefiles
# (i.e., coremake's "make.master"), and are not intended to be used by
# users' makefiles
NCBI_MT_OTHERLIBS = -lthread

NCBI_OPTIONAL_LIBS = BLIB42=libctutils.a BLIB43=libosutils.a

NCBI_THREAD_OBJ = ncbithr.o
NCBI_OTHERLIBS_MT = -lresolv -lsocket -lrpcsvc -lnsl -lgen -lthread -lm
NCBI_THR_OBJ = $(NCBI_LIBDIR)/ncbithr.o
NCBI_THR_ALTOBJ = $(NCBI_ALTLIB)/ncbithr.o
# CodeCenter can't handle the thread library; J. Epstein 5/22/96
NCBI_CLOTHERLIBS = -lresolv -lsocket -lrpcsvc -lnsl -lgen -lm
NCBI_RANLIB = ls -l

NCBI_BIN_MASTER = /netopt/ncbi_tools.master/ver0.0/ncbi/bin
NCBI_BIN_COPY = /netopt/ncbi_tools.copy/ver0.0/ncbi/bin

#NCBI_SYBASE = /netopt/sybase
##NCBI_SYBASE = /netopt/Sybase/clients/10.0.3
##NCBI_SYBLIBS = -L$(NCBI_SYBASE)/lib -R$(NCBI_SYBASE)/lib -lsybdb
##NCBI_SYBFLAG = -I$(NCBI_SYBASE)/include -L$(NCBI_SYBASE)/lib
#NCBI_SYBLIBS = $(NCBI_SYBASE)/libsybdb.a
#NCBI_SYBFLAG = -I$(NCBI_SYBASE)/include 
NCBI_SYBASE = /netopt/Sybase/clients/11.1.0
NCBI_SYBLIBS = -L$(NCBI_SYBASE)/lib  -R$(NCBI_SYBASE)/lib -lsybdb
NCBI_SYBFLAG = -I$(NCBI_SYBASE)/include
NCBI_SYBLIBS_STATIC = $(NCBI_SYBASE)/lib/libsybdb.a

NCBI_SYBLIBS_CT = -L$(NCBI_SYBASE)/lib  -R$(NCBI_SYBASE)/lib -lblk -lct -lcs -ltcl -lcomn -lintl -ltli
#reentrant version:
NCBI_SYBLIBS_CT_r = -L$(NCBI_SYBASE)/lib  -R$(NCBI_SYBASE)/lib -lblk_r -lct_r -lcs_r -ltcl_r -lcomn_r -lintl_r -ltli_r -lthread -ldl

#to compile an Open Server
NCBI_SYBLIBS_OS = -L$(NCBI_SYBASE)/lib  -R$(NCBI_SYBASE)/lib -lsrv -lblk -lct -lcs -ltcl -lcomn -lintl -ltli
#reentrant version:
NCBI_SYBLIBS_OS_r = -L$(NCBI_SYBASE)/lib  -R$(NCBI_SYBASE)/lib -lsrv_r -lblk_r -lct_r -lcs_r -ltcl_r -lcomn_r -lintl_r -ltli_r -lthread -ldl

#NCBI_LAGVIBFLAG = -I/usr/openwin/include -I/netopt/SUNWmotif/include -L/usr/openwin/lib -L/netopt/SUNWmotif/lib -DWIN_MOTIF
NCBI_LAGOTHERLIBS = $(NCBI_OTHERLIBS)
NCBI_LAGVIBFLAG = -I/usr/openwin/include -I/usr/dt/include -L/usr/openwin/lib -L/usr/dt/lib -DWIN_MOTIF
NCBI_DBUGEXTRA = -xsb

#
#FastCGI library for Sun C compilers ver 4.2 and ver 5.0
LIBFASTCGI=-lfcgi`sh -c 'CC -V 2>&1'|cut -f4 -d' '`
