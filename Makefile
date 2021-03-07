#
#  Makefile for mdh, the high energy molecular dynamics code
#  the reed version, and its subprograms (mdsetup, zbl96)
#
#  BEFORE YOU ATTEMPT COMPILING, EDIT the file "config" to contain your
#  machine-dependent settings. Then just say make and hope for
#  the best.
#

#include config

all:	
	cd mdh; make; cd ..
	cd mdhreed; make; cd ..
	./writeheader
	cd mdsetup; make; cd ..
	cd zbl96; make; cd ..

clean:
	- rm */*.o
	- rm */local.h
	- rm zbl96/libzbl96.a
	- rm *~ */*~ */*/~
	- rm core */core */*/core
