
include ../../makefile.in

LIB = ../libpmncbi.a

OBJS =  ncbibase.o \
	ncbilist.o \
	ncbistr.o \
	ncbiio.o

all: $(LIB)

#actually, we don't need ranlib under Solaris, but...
$(LIB): $(OBJS)
	@rm -f $@
	ar q $@ $(OBJS)
	ranlib $@

clean:
	rm -f *.o $(LIB) core *~

install: $(LIB)

purify: $(LIB)

debug: $(LIB)

.KEEP_STATE:

