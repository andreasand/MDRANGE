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

#include "local.h"

#define COLS 54+1
#define ROWS 93+1
#define HEAD 2
#define ACOLS 16
#define BCOLS 38
#define LINE 250

#define max(A,B)  ((A) > (B)) ? (A) : (B)
#define min(A,B)  ((A) < (B)) ? (A) : (B)
