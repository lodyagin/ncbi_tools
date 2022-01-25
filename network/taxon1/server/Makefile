#
#  Makefile for taxon program
#

include $(NCBI)/ncbi.mk
CODEGEN=asncode
ASNTOOL=asntool

PURIFY=/am/purify/purify
CC= $(NCBI_CC) -g

#

INCLUDES=-I$(NCBI_INCDIR)

LIBS1=-lncbitax -lncbiobj -lnetcli -lncbi $(NCBI_OTHERLIBS)
LIBS=-lncbitax1 -lncbiobj -lnetcli -lncbi -ldl $(NCBI_OTHERLIBS) 
LIBDIRS= -L$(NCBI_ALTLIB)

CFLAGS=-g $(INCLUDES)

all : taxserv_fs

# sources and objects for Taxonomy retrieval

taxprocs.o : taxprocs.c
	$(CC) -c $(CFLAGS) -DWARN_ON_FAILURE=1 -o taxprocs.o $(INCLUDES) taxprocs.c

taxsmain.o : taxsmain.c
	$(CC) -c $(CFLAGS) -DWARN_ON_FAILURE=1 -o taxsmain.o $(INCLUDES) taxsmain.c

taxserv_fs: taxprocs.o taxsmain.o
	$(CC) -o taxserv_fs taxsmain.o taxprocs.o $(LIBDIRS) $(LIBS)


asntaxon.l: taxon.asn
	rm -f asntaxon.l*
	$(ASNTOOL) -m taxon.asn -l asntaxon.h
	-mv asntaxon.l* asntaxon.l
	mv asntaxon.h asntaxonstat.h

asntaxon.h: taxon.asn
	$(ASNTOOL) -m taxon.asn -o asntaxon.h

clean:
	rm -f objtaxon.c objtaxon.h asntaxon.l asntaxon.h *.o objtaxonerr objtaxonout








