#
# $Id: hpux_ia64.ncbi.mk,v 1.2 2002/02/19 19:12:32 beloslyu Exp $
#
# testded on:
# HP-UX kiss B.11.20 U ia64 2374655775 unlimited-user license
#
# Contributed by: Balaji Veeraraghavan <balaji_veeraraghavan@hp.com>
#
NCBI_DEFAULT_LCL = hp_ia64
NCBI_MAKE_SHELL = /bin/sh
NCBI_CC = cc -Ae +DD32 +DSitanium -fast -DHPUX -DHPUX_IA64 -Wl,-aarchive_shared
NCBI_CFLAGS1 = -c -z -Wp,-H500000
NCBI_LDFLAGS1 = 
NCBI_OPTFLAG = -fast
NCBI_INCDIR = /usr/ncbi/include
NCBI_LIBDIR = /usr/ncbi/lib
NCBI_ALTLIB = /usr/ncbi/altlib
#will work only when you have Motif installed!
NCBI_VIBFLAG = -I/usr/include/X11R6 -I/usr/dt/include -L/usr/lib/X11R6 -L/usr/dt/lib -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXt -lX11
NCBI_DISTVIBLIBS = -L/usr/lib/X11R6 -L/usr/dt/lib -lXmu -lXm -lXt -lX11
NCBI_OTHERLIBS = -lm
NCBI_RANLIB = ranlib
# Used by makedis.csh
NCBI_MT_OTHERLIBS = -lpthread
NCBI_THREAD_OBJ = ncbithr.o
NETENTREZVERSION = 2.02c2ASN1SPEC6 

NCBI_LBSM_SRC = ncbi_service_lbsmd_stub.c
NCBI_LBSM_OBJ = ncbi_service_lbsmd_stub.o

