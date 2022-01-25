#
# $Id: ibm_pwr3.ncbi.mk,v 1.1 2001/06/15 16:32:10 beloslyu Exp $
#
# Initial version ddas@us.ibm.com 08/19/1999
# Replace r6k with ibm_pwr3 or ibm_auto: cpsosa@us.ibm.com Jun-2001
#
NCBI_DEFAULT_LCL = ibm
NCBI_MAKE_SHELL = /bin/sh
NCBI_CC = xlc_r
#NCBI_CFLAGS1 = -c
NCBI_CFLAGS1 = -c -qcpluscmt
NCBI_LDFLAGS1 = -bmaxdata:0x40000000 -bmaxstack:0x10000000
NCBI_OPTFLAG = -g
NCBI_OPTFLAG = -O3 -qmaxmem=-1 -qarch=pwr3 -qtune=auto -qcache=auto -qcompact -DPOSIX_THREADS_AVAIL
NCBI_INCDIR = /usr/ncbi/include
NCBI_LIBDIR = /usr/ncbi/lib
NCBI_ALTLIB = /usr/ncbi/altlib
#will work only when you have Motif installed!
NCBI_VIBFLAG = -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXmu -lXt -lX11
NCBI_DISTVIBLIBS = -lXm -lXmu -lXt -lX11
NCBI_OTHERLIBS = -lm
#NCBI_OTHERLIBS = -lrs2 -L/usr/local/lib/MASS27 -lmass -lm
NCBI_RANLIB = ranlib
# Used by makedis.csh
NCBI_MT_OTHERLIBS = -lpthread
NCBI_THREAD_OBJ = ncbithr.o
NETENTREZVERSION = 2.02c2ASN1SPEC6 

NCBI_LBSM_SRC = ncbi_service_lbsmd_stub.c
NCBI_LBSM_OBJ = ncbi_service_lbsmd_stub.o
