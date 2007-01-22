CC = gcc-3.4
CPPFLAGS = -I.
CFLAGS = -Wall -g
# -O3 -ffast-math -march=athlon-xp -funroll-loops -fomit-frame-pointer
# -fno-inline-functions
LOADLIBES = -lm
LOADLIBES += -lefence

brontler : brontler.o spr.o init.o io.o utils.o

brontler.o spr.o init.o io.o utils.o: Makefile
brontler.o io.o init.o spr.o: spr.h
io.o init.o spr.o: spr_private.h

.PHONY: clean
clean:
	rm -f *.o brontler
