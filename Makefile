CC = gcc --std=gnu99
CPPFLAGS = -I.
CFLAGS = -Wall -g
# -O3 -ffast-math -march=athlon-xp -funroll-loops -fomit-frame-pointer
# -fno-inline-functions
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
