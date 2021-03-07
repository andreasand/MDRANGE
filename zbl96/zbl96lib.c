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

#include "zbl96lib.h"

double *zbl96(int z1,int z2,double m1,double m2,double rho,
              double min,double max,double step,unsigned int flag)
{
   double scoef[ROWS][COLS];
   double E,x,sunit=1.0,xunit=1.0;
   double *S;
   int i,n;

   readscoef(scoef);
   
   switch(flag & SUNIT){
      case EV_A:
         sunit = 100.0*NA*rho/(m2*1.0e25);      
         break;
      case KEV_NM:
         sunit = NA*rho/(m2*1.0e25);
         break;
      case KEV_UM:
         sunit = 1000.0*NA*rho/(m2*1.0e25);
         break;
      case MEV_MM:
         sunit = 1000.0*NA*rho/(m2*1.0e25);
         break;
      case KEV_UG_CM2:
         sunit = NA/(m2*1e24);
         break;
      case MEV_MG_CM2:
         sunit = NA/(m2*1e24);
         break;
      case KEV_MG_CM2:
         sunit = 1000.0*NA/(m2*1e24);
         break;
      case EV_1E15ATOMS_CM2:
         sunit = 1.0;
         break;
   }
   switch(flag & XUNIT){
      case EV:
         xunit = 1000.0;
         break;
      case KEV:
         xunit = 1.0;
         break;
      case MEV:
         xunit = 0.001;
         break;
      case V0:
         xunit = 1.0;
         break;
      case BETA:
         xunit = 0.0072974;
         break;
      case M_S:
         xunit = 2187673.0;
         break;
      case CM_S:
         xunit = 218767300.0;
         break;
   }

   for(x=min,i=0;x<=max;x+=step,i++);
   n = i;

   S = (double *) malloc(sizeof(double)*n);
  
   for(x=min,i=0;x<=max;x+=step,i++){
      switch(flag & ENERGY){
         case FALSE:
            E = 25.0*x*x/(xunit*xunit);
            break;            
         default:
            E = x/(xunit*m1);
            break;
      }
      switch(z1){
         case 1:
            S[i] = pstop(z2,E,scoef);
            break;
         case 2:
            S[i] = hestop(z2,E,scoef);
            break;
         default:
            S[i] = histop(z1,z2,E,scoef);
            break;
      }
      switch(flag & NUCLEAR){
         case N_ONLY:
            S[i] = nuclear(z1,z2,m1,m2,E*m1);
            break;
         case N_BOTH:
            S[i] += nuclear(z1,z2,m1,m2,E*m1);
            break;
         case N_NO:
            break;            
      }
      S[i] *= sunit;
   }
   
   return(S);
}
