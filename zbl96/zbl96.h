/*
        Zbl96 is a program for calculating electronic and nuclear
        stopping powers according to the semiempirical model of
        Ziegler, Biersack and Littmark.
        
        This program is based on the version 96 of Srim-code.
        
        Version 0.9 written by K. Arstila 16.10.1996

        DO NOT DISTRIBUTE OUTSIDE THE ACCELERATOR LABORATORY OF THE
        UNIVERSITY OF HELSINKI WITHOUT PERMISSION OF THE AUTHOR

                        Kai.Arstila@Helsinki.FI

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "local.h"

#define COLS 54+1
#define ROWS 93+1
#define HEAD 2
#define ACOLS 16
#define BCOLS 38
#define LINE 250

#define NA 6.022e23

#define TRUE  1
#define FALSE 0 

#define MAI 0 
#define NATURAL -1

#define EV_A                 0x0001
#define KEV_NM               0x0002
#define KEV_UM               0x0003
#define MEV_MM               0x0004
#define KEV_UG_CM2           0x0005
#define MEV_MG_CM2           0x0006
#define KEV_MG_CM2           0x0007
#define EV_1E15ATOMS_CM2     0x0008

#define SUNIT                0x000f

#define EV                   0x0010
#define KEV                  0x0020
#define MEV                  0x0030

#define V0                   0x0100
#define BETA                 0x0200
#define M_S                  0x0300
#define CM_S                 0x0400

#define XUNIT                0x0ff0

#define ENERGY               0x00f0
#define VELOCITY             0x0f00

#define N_ONLY               0x1000
#define N_BOTH               0x2000
#define N_NO                 0x3000

#define NUCLEAR              0xf000

#define DSA                  0xf0000

#define DEFAULT (KEV_NM | V0 | N_NO)


#define max(A,B)  ((A) > (B)) ? (A) : (B)
#define min(A,B)  ((A) < (B)) ? (A) : (B)

static char *err_strings[] = {
    "no error",
    "too few command line parameters",
    "maximum energy smaller than minimum energy",
    "negative energy of velocity",
    "no such ion",
    "no such target",
    "no such isotope",
    "negative or zero step",
    "ion velocity exceeds the velocity of light"
};

void readscoef(double [][COLS]);
void readparms(int,char **,int *,int *,double *,double *,double *,double *,
          double *,double *,unsigned int *,double [][COLS]);
double pstop(int,double,double [][COLS]);
double hestop(int,double,double [][COLS]);
double histop(int,int,double,double [][COLS]);
double nuclear(int,int,double,double,double);
void get_element(char *,int,int *,double *,double [][COLS]);
void usage(void);
void fatal_error(int);
