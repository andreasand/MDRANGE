
#
#  Makefile for mdh, the high energy molecular dynamics code, reed version
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


mdhreed:$(SUBS) $(SMALLSUBS) $(HEADERS)
	$(LINK) -o mdhreed $(SUBS) $(SMALLSUBS) -lm
	mv mdhreed ../bin/

$(SMALLSUBS) $(SUBS): $(HEADERS)

