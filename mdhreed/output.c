/* 
   The below defines are needed for automatic defining of CLK_TCK
   in some operating systems, eg. ultrix.
*/
#define _POSIX_SOURCE
#define _INCLUDE_POSIX_SOURCE

#include "common.h"
#include <sys/times.h>
/*
   This file is to include miscellaneous routines used for outputting data
*/
printcoveac(at,nat,t)
struct atom *at;
int nat;
real t;
{
   int i;

   for (i=0;i<nat;i++) printf("co %g %g %g %g %d\n",at[i].x,at[i].y,at[i].z,t,i);
   for (i=0;i<nat;i++) printf("ve %g %g %g %g %d\n",at[i].vx,at[i].vy,at[i].vz,t,i);
   for (i=0;i<nat;i++) printf("ac %g %g %g %g %d\n",at[i].ax,at[i].ay,at[i].az,t,i);
}   
printco(at,nat,box,t,n)
struct atom *at;
int nat;
struct box *box;
real t;
int n;
{
   int i;
   static int c=0;
   
   if (n>0) { 
      if (t==0.0) c=n;
      c++;
      if (c>=n) c=0;
      else return;
   }
   else c=0;

   for (i=0;i<nat;i++) printf("co %.7g %.7g %.7g %g\n",
			      box->movement.x+at[i].x,
			      box->movement.y+at[i].y,
			      box->movement.z+at[i].z,t);
}   

printcoint(at,atflag,nat,box,t,n)
struct atom *at;
struct atomflag *atflag;
int nat;
struct box *box;
real t;
int n;
{
   int i;
   static int c=0;
   
   if (n>0) { 
      c++;
      if (c>=n) c=0;
      else return;
   }
   else c=0;

   for (i=0;i<nat;i++) if (atflag[i].interact) 
     printf("co %.7g %.7g %.7g %g %d\n",box->movement.x+at[i].x,
	    box->movement.y+at[i].y,box->movement.z+at[i].z,t,i);
}   

printrec(recoil,t,n)
struct recoil *recoil;
real t;
int n;
{
   static int c=0;
   
   if (n>0) { 
      if (t==0.0) c=n;
      c++;
      if (c>=n) c=0;
      else return;
   }
   else c=0;

   printf("rec t %g pos %.8g %.8g %.8g vel %g E %g coll %d\n",t,recoil->sx,recoil->sy,recoil->sz,recoil->vel,recoil->energy,recoil->incollision);
   printf("rec2 vx %g vy %g vz %g acc %g eprev %g\n",recoil->vx,recoil->vy,recoil->vz,recoil->acc,recoil->eprev);
   printf("rec3 n %d number %d attype %d\n",recoil->n,recoil->number,recoil->attype);
   printf("rec4 x %g y %g z %g xp %g yp %g zp %g\n",recoil->x,recoil->y,recoil->z,recoil->xp,recoil->yp,recoil->zp);
   printf("rec5 dx %g dy %g dz %g dr %g drp %g\n\n",recoil->dx,recoil->dy,recoil->dz,recoil->dr,recoil->drp);
   oflush();
}   

eltime()
/*  
   Subroutine prints elapsed time since first call
*/
{
   static time_t t1;
   static logical firsttime=True;
   struct tms timebuf;

   if (firsttime) {
      t1=time(NULL);
   }

#ifndef CLK_TCK
#define CLK_TCK 1     /* If nothing else has defined CLK_TCK, use this */
#endif

#define ELTIME 1

#ifdef SunOS4
#undef ELTIME
#define ELTIME 0
#endif

#if (ELTIME == 1)
   times(&timebuf);
   printf("real run time %g sec CPU time %d\n",
      difftime(time(NULL),t1),timebuf.tms_utime/CLK_TCK);
#else
   if (firsttime) printf(" (times() not available :-( )\n");
   else printf("\n");
#endif

   oflush();
   firsttime=False;
}

oflush()
/* 
  Since the standard fflush(NULL) causes problems on some machines,
  use an own that flushes only the needed buffers
*/
{
   fflush(stdout);
   fflush(stderr);
}





