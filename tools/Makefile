# /distrib/tools
#

include ../top.mk

## make BLAST library
##

blast.o : blast.c blast.h blastdef.h blastpri.h
	$(CC) $(CFLAGS) blast.c

blast.ll : blast.c
	$(LINT) $(CFLAGS) blast.c

blastkar.o : blastkar.c blastkar.h 
	$(CC) $(CFLAGS) blastkar.c

lookup.o : lookup.c lookup.h
	$(CC) $(CFLAGS) lookup.c

readdb.o : readdb.c readdb.h
	$(CC) $(CFLAGS) readdb.c

blastutl.o : blastutl.c blastpri.h blastdef.h blast.h
	$(CC) $(CFLAGS) blastutl.c
#$(CC) $(CFLAGS) blastutl.c

gapxdrop.o : gapxdrop.c gapxdrop.h
	$(CC) $(CFLAGS) gapxdrop.c

posit.o : posit.c posit.h
	$(CC) $(CFLAGS) posit.c

dust.o : dust.c dust.h
	$(CC) $(CFLAGS) dust.c

simutil.o : simutil.c simutil.h
	$(CC) $(CFLAGS) simutil.c

ncbisam.o : ncbisam.c ncbisam.h ncbisami.h
	$(CC) $(CFLAGS) ncbisam.c

ncbisort.o : ncbisort.c ncbisort.h ncbisrti.h
	$(CC) $(CFLAGS) ncbisort.c

libblast.a: blast.o blastkar.o lookup.o readdb.o gapxdrop.o blastutl.o posit.o ncbisam.o ncbisort.o dust.o simutil.o blastool.o
	$(AR) $(ARFLAGS) libblast.a blast.o blastkar.o lookup.o readdb.o gapxdrop.o blastutl.o posit.o ncbisam.o ncbisort.o dust.o simutil.o blastool.o

blast2.o : blast2.c
	$(CC) $(CFLAGS) blast2.c
blstout2.o : blstout2.c
	$(CC) $(CFLAGS) blstout2.c

## make applications
##

fastacmd: fastacmd.c libblast.a commindx.c
	$(CC) -o fastacmd fastacmd.c commindx.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB2) $(LIB1) $(OTHERLIBS)

updateindex: updateindex.c commindx.c commindx.h
	$(CC) -o updateindex updateindex.c commindx.c ncbimem.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

fmerge: fmerge.c
	$(CC) -o fmerge fmerge.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB2) $(LIB1) $(OTHERLIBS)

formatdb: formatdb.c
	$(CC) -o formatdb formatdb.c ncbisam.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

formatdb.pure: formatdb.c readdb.h
	purify $(PURECC) -o formatdb.pure formatdb.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB2) $(LIB1) $(OTHERLIBS)

blastall: blastall.c blastutl.c blast.c
	$(CC) -o blastall blastall.c blastutl.c blast.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)
	cp blastall /home/egorov/bin

blastall.pure: blastall.c
	purify $(CC) -o blastall.pure blastall.c blastutl.c blast.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIBDIR)/ncbithr.o $(LIB23) $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)

blasttest: blasttest.c libblast.a dust.o
	$(CC) -o blasttest blasttest.c libblast.a dust.o $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

blasttest.pure: blasttest.c libblast.a dust.o
	purify $(PURECC) -o blasttest.pure blasttest.c libblast.a dust.o $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

newnbrxx: newnbrxx.c libblast.a
	$(CC) -o newnbrxx newnbrxx.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastnew: blastnew.c libblast.a
	$(CC) -o blastnew blastnew.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

blastnew.pure: blastnew.c libblast.a
	purify $(PURECC) -o blastnew.pure blastnew.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastnew.quantify: blastnew.c 
	quantify $(PURECC) -o blastnew.quantify blastnew.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastold: blastold.c libblast.a
	$(CC) -o blastold blastold.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

blastold.pure: blastold.c libblast.a
	purify $(PURECC) -o blastold.pure blastold.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

blsttest: blsttest.c blastutl.c
	$(CC) -o blsttest blsttest.c blastutl.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastngp: blastngp.c libblast.a
	$(CC) -o blastngp blastngp.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

blastngp.quantify: blastngp.c libblast.a dust.o txalign.c txalign.h txdfline.c txdfline.h
	quantify $(PURECC) -o blastngp.quantify blastngp.c libblast.a txalign.c txdfline.c dust.o $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

blastngp.test: blastngp.test.c libblast.a dust.o txalign.c txalign.h txdfline.c txdfline.h
	$(CC) -o blastngp.test blastngp.test.c libblast.a txalign.c txdfline.c dust.o $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB2) $(LIB1) $(OTHERLIBS)

blastngp.pure: blastngp.c libblast.a
	purify $(PURECC) -o blastngp.pure blastngp.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastpgp: blastpgp.c libblast.a
	$(CC) -o blastpgp blastpgp.c txalign.c posit.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastpgp_test: blastpgp_test.c libblast.a
	$(CC) -o blastpgp_test blastpgp_test.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastpgp.pure: blastpgp.c libblast.a
	purify $(PURECC) -o blastpgp.pure blastpgp.c libblast.a $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastpgp.quantify: blastpgp.c libblast.a
	quantify $(PURECC) -o blastpgp.quantify blastpgp.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(NCBI_ALTLIB)/ncbithr.o libblast.a $(LIB23) $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)

blastall.quantify: blastall.c libblast.a
	quantify $(PURECC) -o blastall.quantify blastall.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB23) $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp: blastp.c libblast.a
	$(CC) -o blastp blastp.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp_test: blastp_test.c libblast.a
	$(CC) -o blastp_test blastp_test.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp_test.quantify: blastp_test.c libblast.a
	quantify $(PURECC) -o blastp_test.quantify blastp_test.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)

blp2seqal: blp2seqal.c libblast.a
	$(CC) -o blp2seqal blp2seqal.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB18) $(LIB2) $(LIB1) $(OTHERLIBS)

readdb2fasta: readdb2fasta.c readdb.c readdb.h
	$(CC) -o readdb2fasta readdb2fasta.c readdb.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

flyblast: flyblast.c blastutl.c blast.c blastkar.c
	$(CC) -o flyblast flyblast.c blastutl.c blast.c blastkar.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

flyblast.pure: flyblast.c blastutl.c blast.c blastkar.c
	purify $(PURECC) -o flyblast.pure flyblast.c blastutl.c blast.c blastkar.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

bl2seq: bl2seq.c blastutl.c gapxdrop.c
	$(CC) -o bl2seq bl2seq.c blastutl.c gapxdrop.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

bl2seq.pure: bl2seq.c blastutl.c tofasta.c
	purify $(CC) -o bl2seq.pure bl2seq.c blastutl.c tofasta.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

bl2seq.quantify: bl2seq.c libblast.a
	quantify $(CC) -o bl2seq.quantify bl2seq.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp.pure: blastp.c libblast.a
	purify $(PURECC) -o blastp.pure blastp.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB21) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp.quantify: blastp.c libblast.a
	quantify $(PURECC) -o blastp.quantify blastp.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB21) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp1: blastp1.c libblast.a
	$(CC) -o blastp1 blastp1.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB21) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp1.pure: blastp1.c libblast.a
	purify $(PURECC) -o blastp1.pure blastp1.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB21) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp1.quantify: blastp1.c libblast.a
	quantify $(PURECC) -o blastp1.quantify blastp1.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB21) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp_cmp: blastp_cmp.c libblast.a
	$(CC) -o blastp_cmp blastp_cmp.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB21) $(LIB2) $(LIB1) $(OTHERLIBS)

blastp_cmp.pure: blastp_cmp.c libblast.a
	purify $(PURECC) -o blastp_cmp.pure blastp_cmp.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include libblast.a $(LIB21) $(LIB2) $(LIB1) $(OTHERLIBS)

bldmatrix.o: bldmatrix.c dotmatrx.c dotmatrx.h libblast.a
	$(CC) $(CFLAGS_VIB) bldmatrix.c

bldmatrix:	bldmatrix.o
	$(CC) -o bldmatrix $(LDFLAGS_VIB) bldmatrix.o dotmatrx.c libblast.a $(LIB20) $(LIB23) $(LIB4) $(LIB5) $(LIB18) $(LIB6) $(LIB3) $(LIB2) $(LIB1) $(VIBLIBS) $(OTHERLIBS)

bldmatrix.pure:	bldmatrix.o
	purify $(PURECC) -o bldmatrix.pure $(LDFLAGS_VIB) bldmatrix.o dotmatrx.o libblast.a $(LIB23) $(LIB4) $(LIB5) $(LIB18) $(LIB6) $(LIB3) $(LIB2) $(LIB1) $(VIBLIBS) $(OTHERLIBS)

dotmatm.o : dotmatm.c dotmatrx.h
	$(CC) $(CFLAGS_VIB) dotmatm.c

dotmatrx.o : dotmatrx.c dotmatrx.h
	$(CC) $(CFLAGS_VIB) dotmatrx.c

bldtmatrix.o : bldtmatrix.c
	$(CC) $(CFLAGS_VIB) bldtmatrix.c

dotmatrx:	dotmatm.o dotmatrx.o
	$(CC) -o dotmatrx $(LDFLAGS_VIB) dotmatm.o dotmatrx.o $(LIB23) $(LIB4) $(LIB5) $(LIB18) $(LIB6) $(LIB3) $(LIB2) $(LIB1) $(VIBLIBS) $(OTHERLIBS)

dotmatrx.pure:	dotmatm.o dotmatrx.o
	purify $(PURECC) -o dotmatrx.pure $(LDFLAGS_VIB) dotmatm.o dotmatrx.o $(LIB23) $(LIB4) $(LIB5) $(LIB18) $(LIB6) $(LIB3) $(LIB2) $(LIB1) $(VIBLIBS) $(OTHERLIBS)

bldtmatrix:	bldtmatrix.o libblast.a
	$(CC) -o bldtmatrix $(LDFLAGS_VIB) bldtmatrix.o libblast.a $(LIB23) $(LIB4) $(LIB5) $(LIB18) $(LIB6) $(LIB3) $(LIB2) $(LIB1) $(VIBLIBS) $(OTHERLIBS)

rand_seq: rand_seq.c blastkar.h blastkar.c
	$(CC) -o rand_seq rand_seq.c blastkar.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB2) $(LIB1) $(OTHERLIBS)

gaptest: gaptest.c gapalign.c end2.c g_any.c
	$(CC) -o gaptest gaptest.c gapalign.c end2.c g_any.c $(LDFLAGS) -I. -I/usr/openwin/include -I/netopt/SUNWmotif/include $(LIB23) $(LIB2) $(LIB1) $(OTHERLIBS)

clean	:
	/bin/rm -f  *.o core

backup	:
	/bin/mv BACKUP TAR
	/bin/mv tar.togb.old tar.togb.oldest
	/bin/mv tar.togb tar.togb.old
	tar cvf tar.togb TAR
	/bin/mv TAR BACKUP
	/bin/cp *.prt $(SRC) $(INC) $(SRCCLI) $(INCCLI) Make* BACKUP
