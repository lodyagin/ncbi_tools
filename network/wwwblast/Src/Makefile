# Set the NCBI variable to your local path to the NCBI toolkit!
#NCBI = 

NCBI_INCDIR = $(NCBI)/include
NCBI_ALTLIB = $(NCBI)/altlib
NCBI_LIBDIR = $(NCBI)/lib

NCBI_OTHERLIBS = -lm
THREAD_LIBS = -lpthread

# Uncomment the line(s) relevant to your platform below:

# Solaris Intel 
#NCBI_OTHERLIBS += -lposix4 -lresolv -lsocket -lrpcsvc -lnsl -lgen
# Solaris Sparc
#NCBI_OTHERLIBS += -lposix4 -lresolv -lsocket -lrpcsvc -lnsl
# Darwin
#NCBI_OTHERLIBS = -lc
#THREAD_LIBS =
#DEFINE = -DDARWIN

#
CC=cc 

BLAST_LIBS = -lncbitool -lncbiobj -lncbi $(NCBI_OTHERLIBS)
# For standalone model
ALL_LIBS = -lncbidesk -lvibgif $(BLAST_LIBS) $(THREAD_LIBS)
# For NCBI Client-server model
ALL_LIBS_CS = -lncbitxc2 -lncbiid1 -lnetblast -lnetcli $(ALL_LIBS)
BL2_LIBS = $(BLAST_LIBS)
# For Blast 2 sequences client/server
BL2_LIBS_CS = -lncbiid1 -lnetcli $(BLAST_LIBS)

BL2_DEFINE = -DBL2SEQ_STANDALONE

# For NCBI compilation with signal handling
# The following 2 lines are for internal NCBI BLAST 2 sequences only.
# Comment out for standalone.
#BL2_OBJ = qblastnet.o objBlobj.o $(LIBS)/ncbithr.o
#BL2_LIBS_CS = -lncbiid1 -lnetblast -lnetcli $(BLAST_LIBS) $(THREAD_LIBS)

#DEBUG compilation
LIBS = $(NCBI_ALTLIB)
DEBUG_FLAG = -g

#OPTIMIZED compilation
LIBS = $(NCBI_LIBDIR)
DEBUG_FLAG = -O

OBJ_FILES = wwwbutl.c $(LIBS)/ncbithr.o

#  Defines: 
#  NCBI_CLIENT_SERVER - full NCBI Client/server including BLAST search
#  NCBI_ENTREZ_CLIENT - Client server for gi/accession lookups

INCDIR = -I. -I$(NCBI_INCDIR)

# Additional include path for internal NCBI BLAST 2 sequences only.
# Comment out for standalone
#BLASTQUE = /home/dondosha/BlastQueue
#INCDIR = -I$(BLASTQUE)/Client -I$(BLASTQUE)/ASN.1 -I. -I$(NCBI_INCDIR)

# For standalone model
CFLAGS=  $(DEBUG_FLAG) $(INCDIR)

# Use client/server for gi/accession lookups
CFLAGS_CS = $(DEBUG_FLAG) $(INCDIR) -DNCBI_ENTREZ_CLIENT

# For full NCBI Client-server model
#CFLAGS_FCS= -c $(DEBUG_FLAG) $(INCDIR) -DNCBI_CLIENT_SERVER

.c.o:   Makefile
	$(CC) -c $(CFLAGS) $<

all: psiblast.REAL psiblast_cs.REAL  blast.REAL blast_cs.REAL nph-viewgif.cgi wblast2.REAL wblast2_cs.REAL bl2bag.cgi 

blast.REAL: wwwblast.c $(OBJ_FILES) Makefile
	$(CC)  $(CFLAGS) wwwblast.c $(OBJ_FILES) \
	-L$(LIBS) $(ALL_LIBS) -o ./blast.REAL
	mv ./blast.REAL ..
blast_cs.REAL: wwwblast.c $(OBJ_FILES) Makefile
	$(CC) $(CFLAGS_CS) wwwblast.c $(OBJ_FILES) \
	-L$(LIBS) $(ALL_LIBS_CS) -o ./blast_cs.REAL
	mv ./blast_cs.REAL ..
psiblast.REAL: psiblast.o $(OBJ_FILES) Makefile
	$(CC) $(CFLAGS) psiblast.o $(OBJ_FILES) \
	-L$(LIBS) $(ALL_LIBS) -o ./psiblast.REAL
	mv ./psiblast.REAL ../psiblast.REAL
psiblast_cs.REAL: psiblast.o $(OBJ_FILES) Makefile
	$(CC) $(CFLAGS_CS) psiblast.o $(OBJ_FILES) \
	-L$(LIBS) $(ALL_LIBS_CS) -o ./psiblast_cs.REAL
	mv ./psiblast_cs.REAL ../psiblast_cs.REAL
nph-viewgif.cgi: viewgif.o
	$(CC) viewgif.o -o nph-viewgif.cgi
	mv nph-viewgif.cgi ..
wblast2.REAL : wblast2.c 
	$(CC) $(CFLAGS) $(BL2_DEFINE) wblast2.c \
	-L$(LIBS) $(BL2_LIBS) -o ./wblast2.REAL
	mv wblast2.REAL ..
wblast2_cs.REAL : wblast2.c $(BL2_OBJ)
	$(CC) $(CFLAGS_CS) $(BL2_DEFINE) wblast2.c $(BL2_OBJ) \
	-L$(LIBS) $(BL2_LIBS_CS) -o ./wblast2_cs.REAL
	mv wblast2_cs.REAL ..
bl2bag.cgi : bl2bag.o
	$(CC) bl2bag.o -L$(LIBS) -lncbi -lm -o bl2bag.cgi
	mv bl2bag.cgi ..

xml:
	cd XML && make

clean: 
	rm *.o


