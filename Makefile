CFLAGS=-g -Wall -Ofast
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz
PROG=parse-data kc-c1-fast

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)

parse-data:parse-data.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

kc-c1-fast:kc-c1-fast.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	rm -fr *.dSYM $(PROG)
