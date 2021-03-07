
#include "common.h"
#include <sys/types.h>
#include <unistd.h>

setfiledef(file)
struct file *file;
{
   file->fpo=stdout; strcpy(file->outf,"stdout");
   file->fpdb=stderr; strcpy(file->debugf,"stderr");
   fpdebug=file->fpdb;
   strcpy(file->param,"param.in");
   strcpy(file->coord,"coords.in");
   strcpy(file->elstop,"elstop.in");
   strcpy(file->attrpot,"attrpot.in");
   strcpy(file->reppot,"reppot.in");
   strcpy(file->startf,"startdata.out");
   strcpy(file->rangef,"range.out");
   strcpy(file->depenf,"depen.out");
   strcpy(file->recspecf,"recspec");
   strcpy(file->recspec_Thetaf,"recspec_angle"); /*(JARKKO)*/
   strcpy(file->stopf,"stop.out");
   strcpy(file->avvelf,"avvel.out");

   strcpy(file->cdens,"cdens.in");
   strcpy(file->crho,"charge.in");
   strcpy(file->stopp,"stopping.in");

   strcpy(file->estatef,"estate.out");
   strcpy(file->estatff,"estatfii.out");
   strcpy(file->estattf,"estattheta.out");
}
/*--------------------------------------------------------------------------*/
handleargs(argc,argv,file,reccalc)   /* Handle command line arguments */
char *argv[];
int argc;
struct file *file;
struct reccalc *reccalc;
{
   int i;
   setfiledef(file);

   file->printrecend=False;
   file->printcoordsi=False;
   file->printcoords=False;
   file->printrec=False;
   file->printdepen=False;
   file->printstop=False;
   file->printavvel=False;
   file->fullcalc=False;
   file->printdens=False;
   reccalc->Estat=False;
   reccalc->printz=100000000; /*default large enough*/

   for (i=1;i<argc;i++) {
      if (!strncmp(argv[i],"-h",2)) {
         fprintf(stderr,"Arguments: \n");
         fprintf(stderr,"     -h           help\n");
         fprintf(stderr,"     -d           debug on\n");
         fprintf(stderr,"     -f           always do full calculation of interactions\n");
         fprintf(stderr,"     -nf          don't do full calc\n");
         fprintf(stderr,"     -e           print end result of each recoil\n");
         fprintf(stderr,"     -FD          calculate deposited energies\n");
         fprintf(stderr,"     -nFD         don't calc. deposited energies (default)\n");
         fprintf(stderr,"     -recspec     calculate primary recoil spectrum\n");
         fprintf(stderr,"     -nrecspec    don't calculate primary recoil spectrum (default)\n");
         fprintf(stderr,"     -S           calculate stopping tables \n");
         fprintf(stderr,"     -nS          don't calc. stopping tables (default)\n");
         fprintf(stderr,"     -AV          calculate average velocity (t) tables \n");
         fprintf(stderr,"     -nAV         don't calc. stopping average velocity (t) (default)\n");
         fprintf(stderr,"     -Es          do statistics of E(z) etc.\n");
         fprintf(stderr,"     -nEs         don't do statistics of E(z) etc.\n");
         fprintf(stderr,"     -r           print recoil coordinate, vel and energy\n");       
         fprintf(stderr,"     -p           print coordinates, velocities and accelerations to stdout\n");         
         fprintf(stderr,"     -pi          print coordinates, velocities and accelerations to stdout in inistate calc.\n");         
	 fprintf(stderr,"     -nr          don't  print recoil coordinates\n");
	 fprintf(stderr,"     -np          don't  print coordinates\n");
	 fprintf(stderr,"     -npi         don't  print coordinates in ini calc.\n");
         fprintf(stderr,"     -dm          maximum debug on\n");
         fprintf(stderr,"     -o outfile   redirect stdout\n");
         fprintf(stderr,"     -err errfile redirect stderr\n");
         fprintf(stderr,"     -eo file     redirect stdout and stderr\n");
         fprintf(stderr,"     -pa file     select parametres file\n");
         fprintf(stderr,"     -rp file     select repulsive potential file\n");
         fprintf(stderr,"     -FDf name    change depen file name\n");
         fprintf(stderr,"     -rsf name    change recspec file name\n");
	 fprintf(stderr,"     -rselim min max option to -recscec, energyrange(default=0,E0)\n");
         fprintf(stderr,"     -Sf name     change stopping tables filename\n");
	 fprintf(stderr,"     -pr z        print simulation box coordinates at depth z\n"); /*JARKKO*/
	 fprintf(stderr,"     -dens        print el. dens. spectrum that ion experiences\n");
         exit(0);
      }
      if (!strncmp(argv[i],"-d\0",3)) debug=True;

      if (!strncmp(argv[i],"-f\0",3)) file->fullcalc=True;
      if (!strncmp(argv[i],"-nf\0",4)) file->fullcalc=False;

      if (!strncmp(argv[i],"-e\0",4)) file->printrecend=True;

      if (!strncmp(argv[i],"-FD\0",4)) {
	 file->printdepen=True;
         puts("Output of deposited energies selected");         
      }
      if (!strncmp(argv[i],"-nFD\0",5)) file->printdepen=False;

      if (!strncmp(argv[i],"-recspec\0",4)) {
	 file->recoilspec=True;
         puts("Output of primary recoil spectrum  selected");
         reccalc->ETmin=0;reccalc->ETmax=reccalc->E0;  /*(JARKKO) pohjustus*/
      }
      if (!strncmp(argv[i],"-nrecspec\0",5)) file->recoilspec=False;

      if (!strncmp(argv[i],"-S\0",3)) {
	 file->printstop=True;
         puts("Output of stopping power tables selected");         
      }
      if (!strncmp(argv[i],"-nS\0",4)) file->printstop=False;

      if (!strncmp(argv[i],"-AV\0",4)) {
	 file->printavvel=True;
         puts("Output of average velocities selected");         
      }
      if (!strncmp(argv[i],"-nAV\0",5)) file->printavvel=False;


      if (!strncmp(argv[i],"-Es\0",4)) {
	 reccalc->Estat=True;
         puts("Calculation of recoil E(z) etc. selected");         
      }
      if (!strncmp(argv[i],"-nEs\0",5)) reccalc->Estat=False;

      if (!strncmp(argv[i],"-pi\0",4)) file->printcoordsi=True;
      if (!strncmp(argv[i],"-npi\0",5)) file->printcoordsi=False;

      if (!strncmp(argv[i],"-p\0",3) && strncmp(argv[i],"-pi\0",4)) file->printcoords=True;
      if (!strncmp(argv[i],"-np\0",4) && strncmp(argv[i],"-npi\0",5)) file->printcoords=False;

      if (!strncmp(argv[i],"-r\0",3)) file->printrec=True;
      if (!strncmp(argv[i],"-nr\0",4)) file->printrec=False;
      DEBUG(file->printrec); 

      /*JARKKO*/
      if (!strncmp(argv[i],"-dens\0",3)) file->printdens=True;


      if (!strncmp(argv[i],"-dm\0",4)) {
	 debug=True; maxdebug=True; 
         puts("Maximum debug mode selected !");
      }

      if (!strncmp(argv[i],"-o\0",3)) {
         i++; strncpy(file->outf,argv[i],120);
         file->fpo=fopen(file->outf,"w");
      }
      if (!strncmp(argv[i],"-err\0",5)) {
         i++; strncpy(file->debugf,argv[i],120);
         file->fpdb=fopen(file->debugf,"w");
      }
      if (!strncmp(argv[i],"-eo\0",4)) {
         i++; strncpy(file->debugf,argv[i],120);
         file->fpdb=fopen(file->debugf,"w");
         strncpy(file->outf,argv[i],120);
         file->fpo=fopen(file->outf,"w");
      }

      if (!strncmp(argv[i],"-pa\0",4)) {
         i++; strncpy(file->param,argv[i],120);
      }
      if (!strncmp(argv[i],"-rp\0",4)) {
         i++; strncpy(file->reppot,argv[i],120);
      }

      if (!strncmp(argv[i],"-FDf\0",5)) {
         i++; strncpy(file->depenf,argv[i],120);
      }

      if (!strncmp(argv[i],"-rsf\0",5)) {
         i++; strncpy(file->recspecf,argv[i],120);	
      }

      if(!strncmp(argv[i],"-rselim\0",8)) { /*(JARKKO)*/
	i++; sscanf(argv[i],"%i",&reccalc->ETmin);
	i++; sscanf(argv[i],"%i",&reccalc->ETmax);
	/*printf("%i\n",reccalc->Elimit);*/
      }

      if(!strncmp(argv[i],"-pr\0",8)) {  /*JARKKO*/
	i++; sscanf(argv[i],"%i",&reccalc->printz);
      }


      if (!strncmp(argv[i],"-Sf\0",4)) {
         i++; strncpy(file->stopf,argv[i],120);
      }

   }

   if (debug) fprintf(fpdebug,"FILE NAMES\n%s %s %s\n%s %s %s\n%s %s %s\n%s %s\n",
		      file->outf,file->debugf,file->param,
		      file->coord,file->elstop,file->attrpot,
		      file->startf,file->rangef,file->depenf,
		      file->stopf,file->reppot);

}
/*--------------------------------------------------------------------------*/
start(version,file)
char version[];
struct file file;
{
#ifndef HP9000
   time_t tp;
   pid_t getpid( void );
#endif
   fprintf(file.fpo,"----- %s -----\n\n",version);
 
#ifndef HP9000
   tp=time(NULL);
   fprintf(file.fpo,"Run started %s",asctime(localtime(&tp) ) );
#ifndef NOSHELL
   fprintf(file.fpo,"Pid %d ",(int) getpid()); printhost(file);
#endif
#endif

}

#ifndef NOSHELL

printhost(file)
struct file file;
{
   char *host;

   host=getenv("HOST");

   if (host != NULL) fprintf(file.fpo,"Host %s\n",host);
   else fprintf(file.fpo,"Host <unknown>\n");

}

#endif

