
#
#  Makefile for mdh, the high energy molecular dynamics code
#
include ../config


#
#  Small subroutines that are not likely to change between different versions 
#  of the program.
#
SMALLSUBS = readvar.o readcoords.o random.o time.o solve.o borders.o elstop.o
#
#  SUBS that may be different between different program versions
#
SUBS = mdh.o accelerations.o start.o readin.o initparams.o inistatecalc.o potential.o \
       recoil.o output.o moveslice.o
#
# Headers and other files which all subroutines depend on
#
HEADERS = common.h structs.h arraysizes.h fmath.h global.h Makefile


mdh:	$(SUBS) $(SMALLSUBS) $(HEADERS)
	$(LINK) -o mdh $(SUBS) $(SMALLSUBS) -lm
	mv mdh ../bin/

mdhp:	$(SUBS) $(SMALLSUBS) $(HEADERS)
	$(PROF) -o mdhp $(SUBS) $(SMALLSUBS) -lm
	mv mdhp ../bin	


$(SMALLSUBS) $(SUBS): $(HEADERS)

