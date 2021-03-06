
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

#include "zbl96stop.h"

void readscoef(double scoef[][COLS])
{
   FILE *fp;
   int i,j;
   char buf[LINE];
   
   fp=fopen(F_SCOEFA,"r");

   fgets(buf,LINE,fp);
   fgets(buf,LINE,fp);   

   for(i=1;i<=ROWS;i++)
      for(j=1;j<=ACOLS;j++)
         fscanf(fp,"%lf",&(scoef[i][j]));
         
   fclose(fp);
   
   fp=fopen(F_SCOEFB,"r");

   fgets(buf,LINE,fp);
   fgets(buf,LINE,fp);   

   for(i=1;i<=ROWS;i++)
      for(j=ACOLS+1;j<=ACOLS+BCOLS;j++)
         fscanf(fp,"%lf",&(scoef[i][j]));
         
   fclose(fp);
      
}

double pstop(int z2,double E,double scoef[][COLS])
{
   double sp,x,pe0,pe,sl,sh,ppow;
   
   pe0 = 10.0;

   if(E>1.0e4){
      x = log(E)/E;
      sp = scoef[z2][17] + (scoef[z2][18]*x) + 
           (scoef[z2][19]*x*x) + (scoef[z2][20]/x);
   } else {
      pe = max(pe0,E);
      sl = scoef[z2][9]*pow((double)pe,(double)scoef[z2][10]) +
           scoef[z2][11]*pow((double)pe,(double)scoef[z2][12]);
      sh = scoef[z2][13]/pow((double)pe,(double)scoef[z2][14]) * 
           log((scoef[z2][15]/pe) + scoef[z2][16]*pe);
      sp = sl*sh/(sl + sh);
      if(E <= pe0) {
         ppow = 0.45;
         if(z2 < 7)
            ppow = 0.35;
         sp *= pow((double)E/pe0,(double)ppow);
      } 
   }

   return(sp);
}
double hestop(int z2,double E,double scoef[][COLS])
{
   double se,he0,he,a,b,heh,sp;
   
   he0 = 1.0;

   he = max(he0,E);
   
   b = log(he);
   a = 0.2865 + 0.1266*b - 0.001429*b*b + 0.02402*b*b*b - 
       0.01135*b*b*b*b + 0.001475*b*b*b*b*b;
   a = min(a,30.0);
   heh = 1.0 - exp(-a);

   he = max(he,1.0);
   a = 1.0 + (0.007 + 0.00005*z2)*exp(-pow((double)7.6 - log(he),(double)2));
   heh *= a*a;
   sp = pstop(z2,he,scoef);
   se = sp*heh*4.0;

   if(E <= he0)
      se *= sqrt(E/he0);
   
   return(se);
}     
double histop(int z1,int z2,double E,double scoef[][COLS])
{
   double yrmin,vrmin,v,vfermi,vr,yr,a,q,lambda0,lambda1,zeta0,zeta;
   double sp,se,hipower,eee,eion,l,vfcorr0,vfcorr1,vmin;
   int j;
   
   yrmin = 0.13;
   vrmin = 1.0;
   vfermi = scoef[z2][7];
  
   v = sqrt(E/25.0)/vfermi;
   if(v < 1.0)
      vr = (3.0*vfermi/4.0)*(1.0 + (2.0*v*v/3.0) - (pow((double)v,(double)4)/15.0));
   else
      vr = v*vfermi*(1.0 + 1.0/(5.0*v*v));

   yr = max(vr/pow((double)z1,(double)0.6667),yrmin);
   yr = max(yr,vrmin/pow((double)z1,(double)0.6667));

   a = -0.803*pow((double)yr,(double)0.3) + 1.3167*pow((double)yr,(double)0.6) + 0.38157*yr + 0.008983*yr*yr;

   a = min(a,50.0);
   
   q = 1.0 - exp(-a);

   q = max(q,0.0);
   q = min(q,1.0);

   for(j=22;j<=39 && q>scoef[93][j];j++);
   j--;

   j = max(j,22);
   j = min(j,38);

   lambda0 = scoef[z1][j];
   lambda1 = (q - scoef[93][j])*(scoef[z1][j+1] - scoef[z1][j])/
                                (scoef[93][j+1] - scoef[93][j]);
   l = (lambda0 + lambda1)/pow((double)z1,(double)0.33333);
   zeta0 = q + (1./(2.0*pow((double)vfermi,(double)2)))*(1.0 - q)*
           log(1.0 + pow((double)4.0*l*vfermi/1.919,(double)2));
   a = log(E);
   a = max(a,0.0);

   zeta = zeta0*(1.0 + (1.0/(z1*z1))*(0.08 + 0.0015*z2)*exp(-pow((double)7.6 - a,(double)2)));
   a = vrmin/pow((double)z1,(double)0.6667);
   a = max(a,yrmin);

   if(yr > a){
      sp = pstop(z2,E,scoef);
      se = sp*pow((double)zeta*z1,(double)2);
      eion = min(E,9999.0);
      for(j=41;j<=53 && eion>=scoef[93][j];j++);
      j--;
      j = max(j,41);
      j = min(j,53);

      vfcorr0 =scoef[z2][j];
      vfcorr1 = (eion - scoef[93][j])*(scoef[z2][j+1] - scoef[z2][j])/
                                      (scoef[93][j+1] - scoef[93][j]);
      se *= (vfcorr0 + vfcorr1);
   } else {
      vrmin = max(vrmin,yrmin*pow((double)z1,(double)0.6667));
      a = pow((double)vrmin,(double)2) - 0.8*pow((double)vfermi,(double)2);
      a = max(a,0.0);
      
      vmin = 0.5*(vrmin + sqrt(a));
      eee = 25.0*vmin*vmin;
      sp = pstop(z2,eee,scoef);
      eion = min(eee,9999.0);

      for(j=41;j<=53 && eion>=scoef[93][j];j++);
      j--;
      
      j = max(j,41);
      j = min(j,53);

      vfcorr0 = scoef[z2][j];
      vfcorr1 = (eion - scoef[93][j])*(scoef[z2][j+1] - scoef[z2][j])/
                                      (scoef[93][j+1] - scoef[93][j]);
      sp=sp*(vfcorr0 + vfcorr1);

      hipower = 0.47;
      if(z1 == 3)
         hipower = 0.55;
      else 
         if(z2 < 7)
            hipower = 0.375;
         else
            if(z1<18 && (z2==14 || z2==32))
               hipower = 0.375;
      se = (sp*pow((double)zeta*z1,(double)2))*pow((double)E/eee,(double)hipower);
   }                                         

   return(se);
}
double nuclear(int z1,int z2,double m1,double m2,double E)
{
   double eps,a,sn;
   
   if(E==0.0)
      return(0.0);
   
   eps = 32.53*m2*E/(z1*z2*(m1 + m2)*(pow((double)z1,(double)0.23) + pow((double)z2,(double)0.23)));
   
   if(eps < 30) { 
      a = (0.01321*pow((double)eps,(double)0.21226) + (0.19593*pow((double)eps,(double)0.5)));
      sn = 0.5*log(1.0 + 1.1383*eps)/(eps + a);
   } else {
      sn = log(eps)/(2.0*eps);
   }
   sn *= z1*z2*m1*8.462/((m1 + m2)*(pow((double)z1,(double)0.23) + pow((double)z2,(double)0.23)));
   
   return(sn);
}
