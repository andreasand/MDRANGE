#include "common.h"

inittime(type,unit,time)
struct type *type;
struct unit *unit;
struct time *time;
{
   int i;

   time->kt=time->Kt;
   time->dtreal=time->Dt*1e-15;  /* Recalculated later from recoil->vel and time->kt */
   time->prevdtreal=time->dtreal;
   time->dtmaxreal=time->Dtmax*1e-15;
   for(i=0;i<MAXTYPE;i++) {
      time->dt[i]=time->Dt*1e-15/unit->time[i];
      time->prevdt[i]=time->dt[i];
   }
}
/*---------------------------------------------------------------------------*/
inittsteps(time,unit)
struct time *time;
struct unit *unit;
{
   int i;

   for (i=0;i<MAXTYPE;i++) {
      time->dt[i]=time->dtreal/unit->time[i];
      time->prevdt[i]=time->dt[i];
      time->dtratio[i]=time->dt[i]/time->prevdt[i];
   }
   time->kt=time->Kt;
}
gettsteps(time,unit,new)
struct time *time;
struct unit *unit;
logical new;
{
   int i;

   if (time->dtreal > time->dtmaxreal) time->dtreal=time->dtmaxreal;

   for (i=0;i<MAXTYPE;i++) {
      time->prevdt[i]=time->dt[i];
      time->dt[i]=time->dtreal/unit->time[i];
      time->dtratio[i]=time->dt[i]/time->prevdt[i];
      if (time->dtratio[i] > 1.1001 && !new) {
	 fprintf(fpdebug,"gettsteps: This shouldn't be possible !!!");
	 fprintf(fpdebug,"%d %g %g %g\n",i,time->dt[i],time->prevdt[i],time->dtratio[i]);
      }
   }
   MAXDEBUGSR("Real time step",time->dtreal);
}

calctsteps(time,unit,maxima,recoil,new)
/*
  In calculating the time step, several factors must be taken in account. The
basic relation delta x = v delta t (presented eg. in my gradu) can be very
well used normally. However, when the atom experiences a strong collision, 
this formula alone yields a too big time step. The strong collisions are very
rare, so decreasing the time step all the time would be pretty redundant.

  Therefore some sort of time step dependence on the force the recoil atom
experiences is needed. The force effect on the displacement can be approximated
by calculating 0.5*a*delta t^2. In a very strong collision, the standard criterion
of movement of 0.1 Å per time step is far too big. Let's take as the criterion
in strong collisions instead 0.001 Å, or rather kt/100. Then the accelearation
time step can be calculated with

           delta t = sqrt( 2* kt/100 / a)

 If this time step is smaller than the standard time step calculated from
the velocity, then it should be used. This may cause dramatic decreases in the
length of the time step, but so be it.

  This wasn't very good, this criterion turned out to lead to longer time
steps also at very low energies, which is very unnecessary.

  Let's try another approach: The energy change per time step is

           delta E = a v delta t.

  If we set a maximum available energy change per time step, then we
get a similar criterion for delta t as the ordinary criterion, i.e.

           delta t = delta E_max / (F v)

  A suitable value for  delta E_max might be 300 eV

  This has the advantage that at low velocities and accelerations this
criterion won't be significant. It _is_ significant exactly when
delta E_max/F < Kt.

  For a summary, see K. Nordlund, Comp. Mat. Sci 3 (1995) 448-

*/
struct time *time;
struct unit *unit;
struct maxima *maxima;
struct recoil *recoil;
logical new;
{

   static real vels[VELSTORESIZE];
   static int c=-1,clast=0,nstored=0;
   real vel,velnow,velold;
   real dtrealacc;

   if (new) { c=-1; clast=0; nstored=0; }

   time->prevdtreal=time->dtreal;

   /* Store steps in array vels 

   nstored++; 
   if (nstored>=VELSTORESIZE) { 
      nstored=VELSTORESIZE;
   }
   c++; 
   if (c>=VELSTORESIZE) c=0;
   clast=c+1; 
   if (clast>=VELSTORESIZE) c=0;

   velnow=vels[c]=recoil->vel;
  
   if (nstored >=VELSTORESIZE) velold=vels[clast];
   else velold=velnow;
   vel=MAX(velnow,velold);
   */
   vel=maxima->vel;


/* Calculate time step in SI units according to the bigger velocity */  

   if (vel!=0.0) time->dtreal=(time->kt*unit->length)/
		   (vel*unit->vel[maxima->vel_attype]);

/* Calculate acceleration time step */

   MAXDEBUGSRRRR("calctsteps:",maxima->acc,
		 unit->force,maxima->vel,unit->vel[maxima->vel_attype]); 

   if (maxima->acc==0.0 || maxima->vel==0.0) dtrealacc=1e30;
   else dtrealacc=time->Et*unit->energy/(maxima->acc*unit->force*maxima->vel*
					 unit->vel[maxima->acc_attype]);

   recoil->incollision=False;
   /*
   printf("dt kt %g Et %g v %g a %g\n",time->dtreal,
	  dtrealacc,maxima->vel*unit->vel[maxima->vel_attype],
	  maxima->acc*unit->vel[maxima->acc_attype]);
   */
   if (dtrealacc < time->dtreal) {
      MAXDEBUGSRRRR("calctsteps: selected strong collision time step",
		    dtrealacc,time->dtreal,
		    maxima->vel,maxima->acc);
      recoil->incollision=True;
      time->dtreal=dtrealacc;
   }
   
   /* Prevent too big increases in the time step ! */
   if (!new && time->dtreal > 1.1*time->prevdtreal) time->dtreal=1.1*time->prevdtreal;

   gettsteps(time,unit,new);
}











