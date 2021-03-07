#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include "fmath.h"
/*
  This header does NOT contain common blocks. It should contain very general
stuff that are needed in all subroutine libraries
*/
#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

typedef int logical;
#define True 1
#define False 0
typedef double real; /* This enables simple changing between single and double prec */

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ( (a) > (b) ? (b) : (a) )
#define WITHIN(a,b,c) ( MIN(MAX((a),(b)),(c)) )
#define oabs(a) ( (a) < 0 ? -(a) : (a) ) 

#ifndef abs
#define abs(a) ( (a) < 0 ? -(a) : (a) ) 
#endif

#define sign(a) ( (a) > 0 ? 1 : (a) < 0 ? -1 : 0 )

#include "arraysizes.h"
#include "structs.h"
#include "global.h"

#define iniind(a,b,c,d) ( (a)*MAXATOM*MAXINISTATE*MAXAMORPH+(b)*MAXATOM*MAXINISTATE+(c)*MAXATOM+(d) )
/*layer-am.state-inistate-atom*/
double gaussianrand();
double unirand();



