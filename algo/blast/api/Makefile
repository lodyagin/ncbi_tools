.PHONY: all install

all:
	cd ../lib && make -f Makefile.blast_lib
	make -f Makefile.blast_app

# Solaris only!
CXX = purify /netopt/forte6u2/bin/CC
blast.pure:
	cd ../lib && make -f Makefile.blast_lib
	make APP=$@ CXX="$(CXX)" -f Makefile.blast_app

clean:
	make -f Makefile.blast_app clean
	rm -fr blast.pure

install: blast
	test -d ../bin || mkdir ../bin
	install blast ../bin/
