#include "common.h"
readelstop(file,elstop)
struct file *file;
struct elstop *elstop;
{
   int i,n;
   char buf[170];
   real a,b;
   char *p;

   file->fpel=fopen(file->elstop,"r");
   if (file->fpel == NULL) fprintf(fpdebug,"Can't read elstop file %s\n",file->elstop);

   i=0;
   while(fgets(buf,160,file->fpel) != NULL) {
      if (buf[0]=='%' || buf[0]=='!' || buf[0] == '\n') continue;
      n=sscanf(buf,"%lf %lf",&a,&b);
      if (n<2) break;
      elstop->v[i]=a; elstop->Setable[i]=b;
      i++;
      if (i >= MAXELSTOP) {
         fprintf(fpdebug,"Too many points in elstop array, trying to manage with less\n");
	 fflush(fpdebug);
         break;
      }
   }
   elstop->n=i;
   fprintf(file->fpo,"Read in %d elstop points, max vel %g\n",elstop->n,elstop->v[i-1]);
}

real getelstop(elstop,unit,recoil,v)
struct elstop *elstop;
struct unit *unit;
struct recoil *recoil;
real v;
{
   int i,ihalf,n;
   int iup,idown;
   real vreal,Se;
   real scale;
   real omega,LSSOmega();

   vreal=v*unit->vel[recoil->attype];
   n=elstop->n;
   /* 
      Get two nearest elements with binary search 
   */
   iup=n-1; idown=0;
   while (True) {
      ihalf=idown+(iup-idown)/2;
      if (ihalf==idown) break;
      if (elstop->v[ihalf] < vreal) idown=ihalf;
      else iup=ihalf;
   }

   /* Calculate scaling factor using the formula scale+dscale*E+sqscale*E*E */
   scale=elstop->scale+recoil->energy*(elstop->dscale+elstop->sqscale*recoil->energy);

   /* Get elstop with a simple linear interpolation */
   Se=elstop->Setable[idown]*(elstop->v[iup]-vreal)/(elstop->v[iup]-elstop->v[idown]);
   Se+=elstop->Setable[iup]*(vreal-elstop->v[idown])/(elstop->v[iup]-elstop->v[idown]);
   Se*=scale;

   return(Se);

}

subelstop(at,time,recoil,elstop,unit,type,layer)
struct atom *at;
struct time *time;
struct recoil *recoil;
struct elstop *elstop;
struct unit *unit;
struct type *type;
struct layer layer;  /* Current layer ! */
{
   int i,ihalf,n;
   static int iup=0,idown=1;
   real v,vreal,deltav,Se,Semod;
   real scale;
   real omega,LSSOmega();

   i=recoil->number;
   v=recoil->vel;
   vreal=v*unit->vel[recoil->attype];
   n=elstop->n;

   if (vreal > elstop->v[n-1]) {
     printf("\nElstop array maximum exceeded: vion %g vmaxelstop %g\n\n",
	    vreal,elstop->v[n-1]);
     exit(0);
   }

   if (!(vreal > elstop->v[idown] && vreal < elstop->v[iup])) {
      /* If the previous borders aren't valid, get two nearest 
         elements with binary search 
      */
      iup=n-1; idown=0;
      while (True) {
	 ihalf=idown+(iup-idown)/2;
	 if (ihalf==idown) break;
	 if (elstop->v[ihalf] < vreal) idown=ihalf;
	 else iup=ihalf;
       }
   }
   /* Calculate scaling factor using the formula scale+dscale*E+sqscale*E*E */
   scale=elstop->scale+recoil->energy*(elstop->dscale+elstop->sqscale*recoil->energy);

   /* Get elstop with a simple linear interpolation */
   Se=elstop->Setable[idown]*(elstop->v[iup]-vreal)/(elstop->v[iup]-elstop->v[idown]);
   Se+=elstop->Setable[iup]*(vreal-elstop->v[idown])/(elstop->v[iup]-elstop->v[idown]);
   Se*=scale;

/*  
   Include straggling of electronic stopping if mode is 1
   See A. Luukkainen and M. Hautala, Radiat. Eff. 59 (1982)
*/
   if (elstop->mode==1 && recoil->dr>0.0) {
      /* printf("elstop %g",Se); fflush(NULL);*/
      omega=LSSOmega(type[recoil->attype].Z,layer.Zmean,recoil->vel,
		     layer.ndensity,recoil->dr,unit,recoil->attype);
      Semod=gaussianrand()*omega/recoil->dr;
      MAXDEBUGSRRR("LSS mod",Se,Semod,omega);
      Se+=Semod; 
   }

   deltav=time->dt[recoil->attype]*Se;
   at[i].vx*=1-deltav/v; 
   at[i].vy*=1-deltav/v;
   at[i].vz*=1-deltav/v;

   elstop->Se=Se;

   MAXDEBUGSRRR("elstop.c",elstop->v[idown],elstop->v[iup],Se);

}   

real LSSOmega(Z1,Z2,v,N,deltar,unit,t)
int Z1;
real Z2,v,N,deltar;
struct unit *unit;  
int t;       /* Type of recoil atom */
/*
   Calculate the Omega of the LSS model (Lindhard, Scharff, Schiott, 
   Mat. Fys. Medd. Dan. Vid. Selsk. 33, 14 (1963))
*/
{
   real x,omegaB,L;
   real Bohrvel=2187700.0;   /* SI units */
   real e=1.602189e-19,eps0=8.854188e-12; 
   real pi=3.1415927;
   real vSI;

   omegaB=(Z1*e/(4*pi*eps0)*sqrt(4*pi*Z2*N*deltar/(unit->length*unit->length)));

   vSI=v*unit->vel[t]; 
   x=vSI*vSI/(Bohrvel*Bohrvel*Z2);

   L=1.0;
   if (x < 2.3) {
      L=0.5*(1.36*sqrt(x)-0.016*sqrt(x*x*x));
   }

   MAXDEBUGSRRRR("LSSOmega",v,x,L,omegaB*sqrt(L));

   return(omegaB*sqrt(L));

}






