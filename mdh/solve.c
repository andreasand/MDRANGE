#include "common.h"
/* 
  The algorithms used are the predictor parts of the Smith-Harrison
  algorithm

  It turns out that solve1 and solve2 take up a large fraction of the
  total CPU time in the simulations. Therefore the ugly storing of the
  time constants in arrays, it makes the simulations a lot faster.

*/
solve1(at,time,atflag,box,nat)
struct atom at[];
struct atomflag atflag[];
struct time *time;
struct box *box;
int nat;
{
   int i,t;
   real dt[MAXTYPE],dt2[MAXTYPE],dt6[MAXTYPE],dtr[MAXTYPE];
   real dt26[MAXTYPE],dtr1[MAXTYPE],dtr2[MAXTYPE],dtr3[MAXTYPE];
   real dtr21[MAXTYPE];
   real dtr326[MAXTYPE],dtr26[MAXTYPE],dtr36[MAXTYPE],dtr216[MAXTYPE];

   for(t=0;t<MAXTYPE;t++) {
      dt[t]=time->dt[t];
      dt6[t]=time->dt[t]/6.0;
      dt2[t]=dt[t]*dt[t];
      dt26[t]=dt2[t]/6.0;

      dtr[t]=time->dtratio[t];
      dtr1[t]=1.0+dtr[t];
      dtr2[t]=dtr[t]*dtr[t];
      dtr3[t]=3.0+dtr[t];
      dtr21[t]=dtr2[t]/dtr1[t];

      dtr326[t]=dtr3[t]*dt26[t];
      dtr26[t]=dtr[t]*dt26[t];
      dtr36[t]=dtr3[t]*dt6[t];
      dtr216[t]=dtr21[t]*dt6[t];
   }

   for(i=0;i<nat;i++) {     /* Loop through all atoms */
      t=at[i].type;
      if (atflag[i].interact || atflag[i].recoil) {
         at[i].x+=at[i].vx*dt[t]+dtr326[t]*at[i].ax-dtr26[t]*at[i].a1x;
	 at[i].y+=at[i].vy*dt[t]+dtr326[t]*at[i].ay-dtr26[t]*at[i].a1y;
	 at[i].z+=at[i].vz*dt[t]+dtr326[t]*at[i].az-dtr26[t]*at[i].a1z;

	 at[i].vx+=dtr36[t]*at[i].ax-dtr216[t]*at[i].a1x;
	 at[i].vy+=dtr36[t]*at[i].ay-dtr216[t]*at[i].a1y;
	 at[i].vz+=dtr36[t]*at[i].az-dtr216[t]*at[i].a1z;
      }

      at[i].a2x=at[i].a1x; at[i].a1x=at[i].ax;
      at[i].a2y=at[i].a1y; at[i].a1y=at[i].ay;
      at[i].a2z=at[i].a1z; at[i].a1z=at[i].az;

   }

}


solve2(at,time,atflag,box,nat)
struct atom at[];
struct atomflag atflag[];
struct time *time;
struct box *box;
int nat;
{
   int i,t;
   real dt[MAXTYPE],dt6[MAXTYPE],dtr[MAXTYPE];
   real dtr321[MAXTYPE],dtr3216[MAXTYPE];

   for (t=0;t<MAXTYPE;t++) {
      dt[t]=time->dt[t];
      dt6[t]=dt[t]/6.0;
      dtr[t]=time->dtratio[t];
      dtr321[t]=(3.0+2.0*dtr[t])/(1.0+dtr[t]);
      dtr3216[t]=dtr321[t]*dt6[t];
   }

   for(i=0;i<nat;i++) {     /* Loop through all atoms */
      t=at[i].type;
      if (!atflag[i].interact && !atflag[i].recoil) continue;

      at[i].vx+=dtr3216[t]*at[i].ax;
      at[i].vy+=dtr3216[t]*at[i].ay;
      at[i].vz+=dtr3216[t]*at[i].az;
/* printf("s2 %g %g\n",at[i].ax,at[i].ay,at[i].az); */
/* printf("s2 %g %g %g %g %g %g\n",at[i].vx,at[i].vy,at[i].vz,dtr3216[t],dt[t],dtr[t]); */
   }

}

subtractcmacc(at,nat)
struct atom at[];
int nat;
{
   int i;
   real axtot,aytot,aztot;

   axtot=aytot=aztot=0.0;
   for (i=0;i<nat;i++) {
      axtot+=at[i].ax;
      aytot+=at[i].ay;
      aztot+=at[i].az;
   }
   axtot/=nat;
   aytot/=nat;
   aztot/=nat;
   for (i=0;i<nat;i++) {
      at[i].ax-=axtot;
      at[i].ay-=aytot;
      at[i].az-=aztot;
   }
}

subtractcmvel(at,nat)
struct atom at[];
int nat;
{
   int i;
   real vxtot,vytot,vztot;

   vxtot=vytot=vztot=0.0;
   for (i=0;i<nat;i++) {
      vxtot+=at[i].vx;
      vytot+=at[i].vy;
      vztot+=at[i].vz;
   }
   for (i=0;i<nat;i++) {
      at[i].vx-=vxtot/nat;
      at[i].vy-=vytot/nat;
      at[i].vz-=vztot/nat;
   }
}

setcelltemp(at,atflag,nat,unit,T)
struct atom *at;
struct atomflag *atflag;
struct unit *unit;
int nat;
real T;
{
   int i,n;
   real v,vx,vy,vz,vsq,vnew,vf;

   v=sqrtf(3.0*T/unit->tempconv);
   if (maxdebug) fprintf(fpdebug,"Setting cell temperature to %g K %g vel. units\n",T,v);
   vsq=0.0; n=0;
   vx=vy=vz=0.0;
   for (i=0;i<nat;i++) {
      if (atflag[i].recoil==True) continue;
      n++;
      vx=at[i].vx;
      vy=at[i].vy;
      vz=at[i].vz;
      vsq+=(vx*vx+vy*vy+vz*vz);
   }

   vsq/=n; vnew=sqrt(vsq); 
   if (vnew!=0.0) vf=v/vnew;
   else vf=1.0;
   MAXDEBUGSR("setcelltemp correction factor",vf);
   for (i=0;i<nat;i++) {
      if (atflag[i].recoil==True) continue;
      at[i].vx*=vf;
      at[i].vy*=vf;
      at[i].vz*=vf;
   }
}

getcelltemp(at,atflag,nat,unit,T,ekin)
/*
                                     3                       1  _   2
   Calculate temperature from E    = - N k  T, where E    =  -  >  v
                               kin   2    B           kin    2  ~

*/
struct atom at[];
struct atomflag atflag[];
struct unit *unit;
int nat;
real *T,*ekin;
{
   int i,n;
   real vx,vy,vz;

   *ekin=0.0; n=0;
   vx=vy=vz=0.0;
   for (i=0;i<nat;i++) {
      if (atflag[i].recoil==True) continue;
      n++;
      vx=at[i].vx;
      vy=at[i].vy;
      vz=at[i].vz;
      *ekin+=vx*vx+vy*vy+vz*vz;
   }
   *ekin/=2.0;
   *T= *ekin*2.0*unit->tempconv/(3.0*n);
}





