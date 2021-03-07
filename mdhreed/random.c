#include "common.h"
/* 
   Subroutines to create uniform and gaussian distributed random numbers
   The routines are taken from old MOLDY fortran subroutines (thus the goto).
   The subroutine should first once be called with a seed number, after
   that with the value zero as the argument.

*/
real unirand(seed)
int seed;
{
   static int seedstore;
   static int ii=127773,jj=2147483647;
   static int i1=16807,i2=2836;
   static double p=4.656612875e-10;
   int k1;

   if (seed==0) seed=seedstore;
   if (seed==-1) { fprintf(fpdebug,"unirand.c seed and seedstore for next call %d %d\n",seed,seedstore); return; }

   k1 = seed/ii;
   seed = i1*( seed - k1*ii) - k1 * i2;
   if (seed<0) seed = seed + jj;

   seedstore=seed;
   return ((real)(seed * p));
}
real gaussianrand()
{
   static int iset=0;
   static real gset;

   real v1,v2,r,fac;

   if (iset==0) {
l1:   v1 = 2*unirand(0)-1;
      v2 = 2*unirand(0)-1;
      r = v1*v1+v2*v2;
      if (r>1) goto l1;
      fac = sqrt(-2.*log(r)/r);
      gset = v1*fac;
      iset = 1;
      return (v2*fac);
   }
   else {
      iset = 0;
      return(gset);
    }
}













