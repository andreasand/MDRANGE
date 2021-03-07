#include "common.h" 
initparams(type,unit,time,gen,file)
struct type *type;
struct time *time;
struct unit *unit;
struct gen *gen;
struct file *file;
{
   if (maxdebug) puts("Initialising units");
   initunits(type,unit,file);

   if (maxdebug) puts("Initialising time");
   inittime(type,unit,time);

   unirand(gen->seed);  /* Initialize random number generator */
}
/*---------------------------------------------------------------------------*/
initunits(type,unit,file)
struct type *type;
struct unit *unit;
struct file *file;
{
/*
  The energy "units" are actually conversion factors to SI units.
  To get SI, take the internal unit and multiply with unit.X.
  To get internal units, take SI and divide with unit.X.
  The unit system is a bit complicated, but explained in appendix
  B of Kai Nordlund Pro gradu.
*/
   int i;
   real e=1.6021892e-19;
   real k_B=1.3806662e-23;
   real u=1.6605655e-27;

/* First basic units that don't depend on the time unit */

   fprintf(file->fps,"\n\nUNIT CONVERSION FACTORS (from SI units) \n\n");

   unit->length=1e-10;
   unit->energy=e;
   unit->force=unit->energy/unit->length;
   unit->tempconv=e/k_B; /* Temperature conversion factor eV -> K */

   fprintf(file->fps,"unit length %g m, unit energy  %g J\n",unit->length,unit->energy);
   fprintf(file->fps,"unit force  %g N, temp conv. eV -> K  %g\n",unit->force,unit->tempconv);

/* Then the mass conversion factors. All internal masses are 1 ! */
   fprintf(file->fps,"\nMASS CONVERSION FACTORS, INTERNAL MASSES == 1\n");
   for(i=0;i<MAXTYPE;i++) {
      unit->mass[i]=type[i].m*u;
      fprintf(file->fps,"atom type %d mass %g kg\n",i,unit->mass[i]);
   }

   fprintf(file->fps,"\nMASS AND TIME DEPENDENT INTERNAL UNITS\n");
/* Then the much messier mass and time dependent units */
   for(i=0;i<MAXTYPE;i++) {
      if (type[i].m != 0.0) {
         unit->vel[i]=9822.66/sqrt(type[i].m);
         unit->time[i]=10.1805e-15*sqrt(type[i].m);
         unit->acc[i]=9.6485e17/type[i].m;
         fprintf(file->fps,"atom type %d unit time %g ",i,unit->time[i]);
         fprintf(file->fps,"unit vel %g unit acc. %g\n",unit->vel[i],unit->acc[i]);
      }
   }
   fflush(file->fps);
}
/*---------------------------------------------------------------------------*/
