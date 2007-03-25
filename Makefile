CC = gcc --std=gnu99
CPPFLAGS = -I.
CFLAGS = -Wall -g -O3
# on Core 2 Duo:
# make CC='gcc --std=gnu99 -m32' CFLAGS='-Wall -march=nocona -O3 -ffast-math -funroll-loops'
# on AMD64:
# make CC='gcc --std=gnu99 -m32' CFLAGS='-Wall -march=k8     -O3 -ffast-math -funroll-loops'
# in any case, be sure to use -O3, since it helps much more than any of the other options

# -m32: 32bit pointers take half the space and memory bandwidth of 64bit.
#  ~1.5x speedup for 16 taxa, for brontler -m1 -T5000 -d3 '(tree)'
# A tree topology (e.g. in the duplicate list) is composed entirely of pointers

# -Wall -O3 -ffast-math -march=athlon-xp -funroll-loops -fomit-frame-pointer
# debug: -fno-inline-functions

# with Sun's compiler:
# make CC='c99 -fast -xarch=native' CFLAGS=''

LOADLIBES = -lm
#LOADLIBES += -lefence

.PHONY: all
all: brontler liballspr.a

brontler : brontler.o liballspr.a
LIBOBJS=dupcheck.o spr.o init.o io.o lcg.o utils.o
liballspr.a: $(LIBOBJS)
	ar r $@ $^
#	$(CC) -shared $(CFLAGS) $(LDFLAGS) $(LOADLIBES) -o $@ $^

brontler.o $(LIBOBJS): Makefile
brontler.o: spr.h
$(LIBOBJS): spr.h

.PHONY: clean
clean:
	rm -f *.o brontler liballspr.a
