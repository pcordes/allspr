CC = gcc
CPPFLAGS = -I.
CFLAGS = -Wall -g
# -O3 -ffast-math -march=athlon-xp -funroll-loops -fomit-frame-pointer
# -fno-inline-functions
LOADLIBES = -lm
#LOADLIBES += -lefence

.PHONY: all
all: brontler liballspr.a

brontler : brontler.o -lallspr
liballspr.a: spr.o init.o io.o utils.o
	ar r $@ $^
#	$(CC) -shared $(CFLAGS) $(LDFLAGS) $(LOADLIBES) -o $@ $^

brontler.o spr.o init.o io.o utils.o: Makefile
brontler.o io.o init.o spr.o: spr.h
io.o init.o spr.o: spr_private.h

.PHONY: clean
clean:
	rm -f *.o brontler
