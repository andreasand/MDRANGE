
#
#  Makefile for zbl96
#

include ../config

CFLAG=
# CFLAG= -g  -pg
# CFLAG= -O6 -fstrength-reduce  -fthread-jumps -frisc

LIB= -lm
LFLAG = # -g -pg -static

OBJS= zbl96stop.o zbl96main.o 
LIBOBJS= zbl96stop.o zbl96lib.o

zbl96: $(OBJS) libzbl96.a
	$(CC) $(LFLAG) -o zbl96 $(OBJS) $(LIB)
	mv zbl96 ../bin
	
lib: libzbl96.a $(LIBOBJS)

clean:
	rm  $(OBJS) zbl96

cleanlib:
	rm  $(LIBOBJS) libzbl96.a testlib testlib.o

testlib: testlib.o $(LIBOBJS) libzbl96.a
	$(CC) $(LFLAG) -o testlib testlib.o libzbl96.a $(LIB)

libzbl96.a: $(LIBOBJS)
	ar r libzbl96.a $(LIBOBJS)
	ranlib libzbl96.a

testlib.o: testlib.c testlib.h local.h Makefile
	$(CC) -c $(CFLAG) testlib.c	

zbl96stop.o: zbl96stop.c zbl96stop.h local.h Makefile
	$(CC) -c $(CFLAG) zbl96stop.c

zbl96main.o: zbl96main.c zbl96.h local.h Makefile
	$(CC) -c $(CFLAG) zbl96main.c

zbl96lib.o: zbl96lib.c zbl96lib.h local.h Makefile
	$(CC) -c $(CFLAG) zbl96lib.c
	
