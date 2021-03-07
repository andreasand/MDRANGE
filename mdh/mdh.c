#define MAIN 

#include "common.h" 

#define VERSION "mdh V3.09c Oct 6 2016"

/*
  Version history (more detailed info in file modification_history):

  0.95:    First about working version, inistate Temp still behaved wildly. 
  0.95.1:  potential stored in array, sigmar and straggling calculation added
  0.96:	   Temperature calculations corrected, inistate calc. now correct !
  0.97:    Morse potential added, other minor improvements (=OMI)
  0.98:    Stopway array added, and OMI.
  0.99:    Layer handling added + OMI.
  0.99.1:  Added time->Et and recoil.incollision criteria to algorithms
  1.0:     Added polycryst. calc. and reccalc->type==2.
  1.01:    Deposited energies calculation added.
  1.01a:   Bugs in surface handling for Trange=1 corrected
  1.1:	   Added atom changing routines (type[].prob) for steel sim.
  1.1a:	   Corrected bug which caused tilts in Linux systems and other bugs.
  1.2:     Elstop linear scaling, thorough testing etc.
  1.2a:    Morse potential values for Ta and Nb changed
  1.3:     Added small angle change routines in polycr. calc.
  1.3b:    Change in initial theta selection + ind. repulsive pot.
  1.3c-d:  Corrected large number of stupid bugs
  1.4:     Minor changes in layer handling routines
  1.5:     Added calculation of straggling of electronic stopping
  1.5a-d:  Minor bug corrections
  1.51:    Stopping power calculation added
  1.51a-f: Bug corrections and changes in stopping power calc.
  1.52:    Added average velocity calculation
  1.52b:   Minor changes in avvel stuff
  1.6x:    Temporary version with some of JS:s modifications
  1.7:     Main distribution version
  1.72:    Changed splinereppot to interpolate screening function
  1.74:    Corrected factor sqrt(3) in Debye displacement calculation
  1.8      Made float->double, improved ZBL and spline pot. interpolations
  1.82:    (JARKKO)-indexed modifications, angular recoil spectrum with energy-range
           option (-rselim min max)
  2.00Beta (JARKKO) added damagebuildup+sputtering+minor optimizations
  3.0Z     rewritten moveslice, added debuggin option -pr z, useful for 
           monitoring the surrounding of the ion at depth z.

  3.02     Corrected error in spline derivative in Jarkkos version
  3.03     Added periodic cell shift option
  3.05     Corrected recspec_angle limit handling. Added -collisionsout output option
  3.07     Corrected projected range surface calculation
  3.09     Added features for antiproton calculations
*/
/*
   Different elstops for different layers not yet implemented for puskastopping.
*/
main(argc,argv)
int argc;
char *argv[]; 
{
   struct unit unit;
   static struct atom *at;          /* Atom coordinates */
   static struct atomex *atex;      /* Miscellaneous atom properties */
   static struct atomflag *atflag;  /* flags etc. */
   struct type type[MAXTYPE];        /* Atom type data */
   struct time time;                 /* time units depends on atom mass */
   static struct iniat *iniat;      /* Initial positions */
   struct layers layers;
   struct potcrit potcrit;
   struct pot pot;                /* pot.rep and pot.attr structures */
   struct elstop elstop;          /* Elstop data */
   struct recoil recoil;
   struct maxima maxima;
   struct reccalc reccalc;        /* Recoil calculation type etc. */
   struct recresult recresult;    /* Recoil calculation result data */
   struct file file;              /* File data */
   struct box box;                /* Simulation cell sizes etc. */
   struct physical physical;      /* Physical constants of the simulation cell */
   struct gen gen;                /* General miscellaneous params */
   struct poly poly;              /* Internal polycryst. calc. params */

   int nrec;      /* Number of recoil event currently calculated */
   int ntstep;    /* Number of time step */
   real t;        /* Time */
   int noutput,nat;        /* Current number of atoms */
   int nlayer,ninistate,namorph;   /* Current layer and inistate calc. selected, namorphstate JARKKO-damage */
   int ngbrlistmode;
   logical update;    /* Force update of neighbour list */
   logical firstrecstep;  /* True if first time step for new recoil */
   logical recerror,recbound;

   double itime;   /* Index of time steps */
   int i;       /* way of stopping */

   /*variables for damage-buildup and sputtering   */
   real Ymean=0.0,ddz=5.43,r2,FD,expath=0;
   int j,k;
   logical templog;
   logical issput;
   struct atom *attemp;
   struct atomflag *atflagtemp;
   struct atomex *atextemp;
   point origmovelim,origboxsize;
   point tat;
   FILE *fad,*fmov;
   int frame=0;
   
/* 
  Note that structures are sent as pointers to subroutines a) when the
  structure contents are changed and b) when it is necessary for speed
*/
   fpdebug=stderr;  /* Temporary value, real given in handleargs */
   eltime();

   handleargs(argc,argv,&file,&reccalc);   /* Handle command line arguments */

   start(VERSION,file);   /* Print date, program version, pid etc */

   if (debug) puts("Going into parameter read");
   readin(type,&time,&layers,&potcrit,&pot,&recoil,&reccalc,&box,
          &physical,&file,&gen,&elstop);
   
   if (file.recoilspec) {
     if (reccalc.ETmin<=0.0 && reccalc.ETmax<=0.0) {
       /* not read in, initialize them now */
       reccalc.ETmin=0.0;reccalc.ETmax=reccalc.E0;    
     }
     printf("Angular recoil specctrum calculated in E interval %g - %g\n",reccalc.ETmin,reccalc.ETmax);
   }
   
   if(physical.Damage || physical.sputtering)
     {
     physical.Dosescaling=physical.Dose*(1e-8)*(1e-8)*(5.43)*(5.43)/reccalc.Ncalc;
     printf("DOSE SCALING in damage buildup model=%g\n",physical.Dosescaling);
     }
   /*JARKKO-damage sc=Dose(1/Ongstrom^2)/(nions/surface_area_of_damage_box)*/
   

/*large atom tables for the simulation cells*/
   at= (struct atom *) malloc(MAXATOM*sizeof(struct atom));
   attemp = (struct atom *) malloc(MAXATOM*sizeof(struct atom)); /*JARKKO-sputtering*/
   atex= (struct atomex *) malloc(MAXATOM*sizeof(struct atomex));
   atextemp= (struct atomex *) malloc(MAXATOM*sizeof(struct atomex)); /*JARKKO-sputtering*/
   atflag= (struct atomflag *) malloc(MAXATOM*sizeof(struct atomflag));
   atflagtemp= (struct atomflag *) malloc(MAXATOM*sizeof(struct atomflag)); /*JARKKO-sputtering*/

   iniat= (struct iniat *) malloc(MAXLAYER*MAXAMORPH*MAXINISTATE*MAXATOM*sizeof(struct iniat));
   MAXDEBUGSR("Initialized iniat array to ",MAXLAYER*MAXAMORPH*MAXINISTATE*MAXATOM*sizeof(struct iniat));

   if (debug) puts("Going into parameter initialization");
   initparams(type,&unit,&time,&gen,&file);

   if (debug) puts("Going into atom coordinate reading");
   readcoords(iniat,&layers,&file,type,&box,&unit,&physical);

   /*  for(i=0;i<MAXTYPE; i++)
    for(j=0;j<ORIGMAXLAYER; j++) 
       printf("%d %d %d\n",i,j,file.layeratoms[i][j]);
   */


   if (debug) puts("Going into inistate calculation");
   inistatecalc(iniat,&pot,&potcrit,&time,&unit,type,&file,&box,&physical,&layers);
   fprintf(file.fpo,"Inistate calculation ended\n");  eltime(); 
   noutput=0;

   fprintf(file.fpo,"\nStarting recoil event calculation\n");
   fprintf(file.fpo,  "---------------------------------\n\n");
   /* Then for the actual recoil calculation loop */

   readelstop(&file,&elstop,&type,&layers); /*Nlayers must be provided 
						to readelstop*/

   initrec(&file,&reccalc,&recresult,&poly,&recoil,&unit,&elstop,type);  
	/*Init some rec. calc. vars */

   if (debug) puts("Starting recoil calculation loop");
   oflush();
   
   origmovelim=box.movelim;  /*original movelim:s*/
   origboxsize=box.size;
   
   if(file.movie) fmov=fopen("at.movie","w");

   if (reccalc.collisionsout>=0.0) {
     printf("Outputting collisions.txt file of all recscpec recoils with E > %g\n",reccalc.collisionsout);
     printf("\n WARNING: THIS MAY RESULT INTO A HUGE FILE!! \n\n");
     file.fpcollisions=fopen("collisions.txt","w");
   }


   for (nrec=1;nrec<=reccalc.Ncalc;) {
      firstrecstep=True;
      
      getlayer(&layers,reccalc.Startmax.z,&nlayer,firstrecstep);

      if(physical.scaling)
	{
	  /*if scaling done in the readcoords.c, the movelims must be set right for the next ion.*/
	if(debug) fprintf(fpdebug,"Bscaling is set on and movelim.sizes are scaled for the first layer\n");
	box.movelim=origmovelim;
 	box.size=origboxsize;
	}
 
      /*get damage state*/
      getamorphstate(&namorph,&recresult,reccalc.Startmax.z,box.movelim.z,&physical,&layers,nlayer); 
     
      nat=layers.layer[nlayer].natoms[namorph]; /*number of atoms*/
      DEBUGSRRR("Starting recoil ",nrec,nat,poly.ndomain);
      /* Zero arrays, set open boundaries ... */ 

      startrec(iniat,at,atex,atflag,&layers,&nlayer,&ninistate,&box,&potcrit,&time,
	       &reccalc,&recresult,&poly,&elstop,&namorph);  /*JARKKO-damage  */

      
      /*check if first layer has different dimensions than simulation box=at[], 0.01 is tolerance*/
      if(fabs(layers.layer[nlayer].boxsize[namorph].x-box.size.x)>0.01 ||
	 fabs(layers.layer[nlayer].boxsize[namorph].y-box.size.y)>0.01 ||
	 fabs(layers.layer[nlayer].boxsize[namorph].z-box.size.z)>0.01)
	{
	printf("%f,%f,%f,%f,%f,%f\n",layers.layer[nlayer].boxsize[namorph].x,box.size.x,
	  layers.layer[nlayer].boxsize[namorph].y,box.size.y,layers.layer[nlayer].boxsize[namorph].z,box.size.z);
	j=0;
        printf("nat before simboxcheck: %d\n",nat);
        for(k=0;k<nat;k++)
	   if(at[k].x>=0.0 && at[k].x<box.size.x) 
	   if(at[k].y>=0.0 && at[k].y<box.size.y) 
	   if(at[k].z>=0.0 && at[k].z<box.size.z)   
	     {copyatom(&at[k],&atex[k],&atflag[k],&attemp[j],&atextemp[j],&atflagtemp[j]);j++;}
      
        for(k=0;k<j;k++) 
           copyatom(&attemp[k],&atextemp[k],&atflagtemp[k],&at[k],&atex[k],&atflag[k]);
	nat=j;
	printf("nat after simboxcheck: %d\n",nat);
	}
	

      setcelltemp(at,atflag,nat,&unit,physical.Trec); 

      t=0.0;
      if (file.printcoords) printco(at,nat,&box,0.0,0);
      MAXDEBUGS("Creating new recoil\n");
      createrecatom(&reccalc,&recoil,&maxima,&nat,at,atflag,&poly,nrec);
      
      if(physical.sput_surf>0)  /* surface modifications*/
        erodesurface(at,atex,atflag,attemp,atextemp,atflagtemp,&recoil,&physical,&box,&nat);

      for(j=0; j<ZBOXES; j++)
       recresult.depennucl_boxtemp[j]=0.0;

      /* Select neigbour list mode, see accelerations.c */
      if (!file.fullcalc) ngbrlistmode=1;   
      else ngbrlistmode=0; 

      time.dtreal=time.Dt*1e-15;
      inittsteps(&time,&unit);
      calctsteps(&time,&unit,&maxima,&recoil,True);

      for(i=0;i<layers.Nlayers;i++) /*setting box movement to zero for all layers.*/
	for(k=0;k<layers.layer[i].Namorph; k++)
	  {
	  layers.layer[i].front[k].x=0.0;
          layers.layer[i].front[k].y=0.0;
          layers.layer[i].front[k].z=0.0;
	  }


      oflush();
      update=True; recerror=False; recbound=False;
      itime=0.0;
      while ( (recoil.energy>reccalc.Emin || recoil.incollision) ) {
	itime+=1.0;

	/* As of V3.07, handling of all varieties of sputtering moved to end of time loop */

	/*if -pr option used*/
	if(recoil.sz>reccalc.printz)
	  {
	    fad=fopen("at.xyz","w");
	    fprintf(fad,"%d\n\n",nat);
	    for(i=0;i<nat-1;i++)
	    {	     
	    if(at[i].type==1) fprintf(fad,"Si %lg %lg %lg %d\n",at[i].x,at[i].y,at[i].z,1);
	    else fprintf(fad,"H %lg %lg %lg %d\n",at[i].x,at[i].y,at[i].z,2);
	    }
	  fprintf(fad,"Cu %lg %lg %lg %d\n",at[nat-1].x,at[nat-1].y,at[nat-1].z,0);
	  reccalc.printz=100000000;
	  fclose(fad);
	  }
	
	/* if movie option used*/
	if(nrec==1)
	if(file.movie && (recoil.r_path-expath)>0.5)
	  {	   
	    expath=recoil.r_path; frame++;
	  fprintf(fmov,"%d\n",nat);
	  fprintf(fmov,"Frame number %d %lg fs\n",frame,time.Dt);
	  tat.x=at[nat-1].x-box.size.x/2;
	  tat.y=at[nat-1].y-box.size.y/2;
	  tat.z=at[nat-1].z-box.size.z/2;
	  for(i=0;i<nat-1;i++)
	    {
	    if(at[i].type==1) fprintf(fmov,"Si %lg %lg %lg %d %d\n",at[i].x-box.size.x/2-tat.x,at[i].y-box.size.y/2-tat.y,
				      at[i].z-box.size.z/2-tat.z,1,i);
	    else fprintf(fmov,"H %lg %lg %lg %d %d\n",at[i].x-box.size.x/2-tat.x,at[i].y-box.size.y/2-tat.y,
			 at[i].z-box.size.z/2-tat.z,2,i);
	    }
	  fprintf(fmov,"Cu %lg %lg %lg %d %d\n",at[nat-1].x-box.size.x/2-tat.x,at[nat-1].y-box.size.y/2-tat.y,
		  at[nat-1].z-box.size.z/2-tat.z,0,nat-1);	  
	  }

	 MAXDEBUG(t);
	 MAXDEBUGS("\nSolving 1\n");
         solve1(at,&time,atflag,&box,nat);
	 MAXDEBUGS("Going into slice move\n");
	 
	 moveslice(at,atex,atflag,type,iniat,&layers,&nlayer,&nat,&box,
		     &recoil,&reccalc,&recresult,&update,firstrecstep,
		     file.recoilspec,&physical,&unit,&file);  /*JARKKO-damage*/

	 MAXDEBUGS("Getting accs\n");
         accelerations(at,atex,atflag,&pot,&potcrit,&box,
		       &physical,&unit,nat,ngbrlistmode,&update,&file);
	 MAXDEBUG(atex[recoil.number].pot);
	 MAXDEBUGS("Solving 2\n");
         solve2(at,&time,atflag,&box,nat);
	 MAXDEBUGS("Getting rec prop\n");
	 getrecoilprop(&recoil,&reccalc,at,&box,&poly); 
	 MAXDEBUGS("Getting maximum velocities etc.\n");
	 getmaximumprop(&maxima,at,nat); 
	 MAXDEBUGS("Subtracting elstop\n");

	 /*printf("Fde:%lg, Fdn:%lg\n",elstop.Se*recoil.dr,(recoil.eprev-recoil.energy)-elstop.Se*recoil.dr);*/
	 /*if(recoil.sz>10000) printf("E: %lg ,z= %lg\n",recoil.energy,recoil.sz);*/
	 if(recoil.sz>=physical.sput_surf) /*if recoil below real surface */
           subelstop(at,&time,&recoil,&elstop,&unit,&type,layers.layer[nlayer],&potcrit,nlayer,&file); 
	 else{elstop.totalfirsov=0;
	      elstop.Se=0;}

  
	 if (reccalc.Estat) getestat(&recoil,&reccalc,&file,False);

/*for damage model, deposited energy calculated continiously*/
	 check_depencalc(&physical,&elstop,&recoil,&reccalc,&box,&poly,&recresult,&layers,0,&Ymean);
	
	 if (file.printdepen || file.printstop || file.printavvel) 
	    calcmisc(&recoil,&elstop,&reccalc,&recresult,&poly,&file,
		      firstrecstep,t/unit.time[recoil.attype],&physical);  
	                   /* ,layers.layer[nlayer].ndensity); */
	 MAXDEBUGS("Getting new t step\n");
         calctsteps(&time,&unit,&maxima,&recoil,False);
	 
	 if (reccalc.poly && recoil.r_proj > poly.cursize) {
	    exitpolydomain(&poly,&recoil);
	    break;   /* Exit loop over time steps */
	 }

	 if (recoil.energy > 1.4* reccalc.E0) {
	   if (reccalc.allowhugeenergy==0) { 
	     fprintf(fpdebug,"Warning: Recoil n. %d has insane energy %g.\n",nrec,recoil.energy);
	     fprintf(fpdebug,"Discarding it... try smaller timestep!\n"); 
	     printrec(&recoil,t,0);
	     printcoint(at,atflag,nat,&box,t,0);
	     recerror=True;   /* Discard this recoil event */
	     break;   /* Exit loop over time steps */
	   }
	   else {
	     if (reccalc.allowhugeenergy==1) { 
	       fprintf(fpdebug,"Warning: Recoil n. %d has insane energy %g.\n",nrec,recoil.energy);
	        fprintf(fpdebug,"Since reccalc.allowhugeenergy=%d continuing anyway\n",reccalc.allowhugeenergy); 
	     }
	     else  if (reccalc.allowhugeenergy==2) { 
	       if (recoil.energy >100* reccalc.E0) {
		 fprintf(fpdebug,"Warning: Recoil n. %d has extremely high energy, 100xinitial %g.\n",nrec,recoil.energy);
		 fprintf(fpdebug,"Discarding it assuming nuclear reaction occurred\n"); 
		 fprintf(file.fpo,"Bound recoil end pos %.2f %.2f %.2f E %.8g\n",recoil.sx,recoil.sy,recoil.sz,recoil.energy);
		 recbound=True;   /* Discard this recoil event */
		 break;   /* Exit loop over time steps */
	       }
	     }
	   }
	 }
	 if (reccalc.allowhugeenergy==3) { 
	   if (physical.rmin < 1e-5) {
	     fprintf(fpdebug,"Warning: Recoil n. %d is closer than 1 fm to atom %g\n",nrec,physical.rmin);
	     fprintf(fpdebug,"Discarding it assuming nuclear reaction occurred\n"); 
	     fprintf(file.fpo,"Bound recoil end pos %.2f %.2f %.2f E %.8g\n",recoil.sx,recoil.sy,recoil.sz,recoil.energy);
	     recbound=True;   /* Discard this recoil event */
	     break;   /* Exit loop over time steps */
	   }
	   if (itime > 1e9) { 
	     fprintf(fpdebug,"Warning: Recoil n. %d has been simulated more than 1 billion steps. Giving up and exiting\n",nrec);
	     fprintf(fpdebug,"Energies: Epot %g Ekin %g\n",atex[recoil.number].pot,physical.ekin);
	     fprintf(file.fpo,"Bound recoil end pos %.2f %.2f %.2f E %.8g\n",recoil.sx,recoil.sy,recoil.sz,recoil.energy);
	     printrec(&recoil,t,0);
	     recbound=True;   /* Discard this recoil event */
	     break;   /* Exit loop over time steps */
	   }	      

	 }
	 if (reccalc.allowhugeenergy==4) { 
	   if (physical.rmin < 1e-4) {
	     fprintf(fpdebug,"Warning: Recoil n. %d is closer than 10 fm to atom %g\n",nrec,physical.rmin);
	     fprintf(fpdebug,"Discarding it assuming nuclear reaction occurred\n"); 
	     fprintf(fpdebug,"Energies: Epot %g Ekin %g\n",atex[recoil.number].pot,recoil.energy);
	     fprintf(file.fpo,"Bound recoil end pos %.2f %.2f %.2f E %.8g\n",recoil.sx,recoil.sy,recoil.sz,recoil.energy);
	     recbound=True;   /* Discard this recoil event */	     
	     break;   /* Exit loop over time steps */
	   }
	   if (itime > 1e8) { 
	     fprintf(fpdebug,"Warning: Recoil n. %d has been simulated more than 100 million steps. Giving up and exiting\n",nrec);
	     fprintf(fpdebug,"Energies: Epot %g Ekin %g\n",atex[recoil.number].pot,recoil.energy);
	     printrec(&recoil,t,0);
	     recbound=True;   /* Discard this recoil event */
	     break;   /* Exit loop over time steps */
	   }	      
	 }
	 if (reccalc.allowhugeenergy==5) { 

	   if (atex[recoil.number].pot < -recoil.energy  ) {
	     fprintf(fpdebug,"Warning: Recoil potential energy %g less than -kinetic %g\n",atex[recoil.number].pot,recoil.energy);
	     fprintf(fpdebug,"Assuming atoms ended up in bound orbit. "); 
	     fprintf(fpdebug,"Interatomic minimum separation %g\n",physical.rmin);
	     fprintf(file.fpo,"Bound recoil end pos %.2f %.2f %.2f E %.8g\n",recoil.sx,recoil.sy,recoil.sz,recoil.energy);
	     recbound=True;   /* Discard this recoil event */	     
	     break;   /* Exit loop over time steps */
	   }
	   if (itime > 1e8) { 
	     fprintf(fpdebug,"Warning: Recoil n. %d has been simulated more than 100 million steps. Giving up and exiting\n",nrec);
	     fprintf(fpdebug,"Energies: Epot %g Ekin %g\n",atex[recoil.number].pot,recoil.energy);
	     printrec(&recoil,t,0);
	     recbound=True;   /* Discard this recoil event */
	     break;   /* Exit loop over time steps */
	   }	      
	 }
	 if (reccalc.allowhugeenergy==6) { 
	   if (physical.rmin < 3e-5) {
	     fprintf(fpdebug,"Warning: Recoil n. %d is closer than 3 fm to atom %g\n",nrec,physical.rmin);
	     fprintf(fpdebug,"Discarding it assuming nuclear reaction occurred\n"); 
	     fprintf(fpdebug,"Energies: Epot %g Ekin %g\n",atex[recoil.number].pot,recoil.energy);
	     fprintf(file.fpo,"Bound recoil end pos %.2f %.2f %.2f E %.8g\n",recoil.sx,recoil.sy,recoil.sz,recoil.energy);
	     recbound=True;   /* Discard this recoil event */	     
	     break;   /* Exit loop over time steps */
	   }
	   if (itime > 1e8) { 
	     fprintf(fpdebug,"Warning: Recoil n. %d has been simulated more than 100 million steps. Giving up and exiting\n",nrec);
	     fprintf(fpdebug,"Energies: Epot %g Ekin %g\n",atex[recoil.number].pot,recoil.energy);
	     printrec(&recoil,t,0);
	     recbound=True;   /* Discard this recoil event */
	     break;   /* Exit loop over time steps */
	   }	      
	 }


	 if (recoil.range>reccalc.Zmax || recoil.range>layers.layer[layers.Nlayers-1].maxz) break;

	 issput=False;
	 if (recoil.sz<(2.0*reccalc.Startmin.z)) { issput=True; break; }
	 /* V3.07: revised these to handle projected range properly */
	 if (reccalc.type!=2) {
	   /* Analyze cutoff criteria by z */
	   if (reccalc.Trange == 1) {
	     if (recoil.range<(2.0*recoil.z0/cos(recoil.theta) )) { issput=True; break; }
	   }
	   else {
	     if (recoil.range<(2.0*reccalc.Startmin.z)) { issput=True; break; }
	   }
	 }

         if (file.printcoords) printcoint(at,atflag,nat,&box,t,reccalc.Outputstep);
         if (file.printrec || maxdebug) printrec(&recoil,t,reccalc.Outputstep);
         t+=time.dtreal;
	 if (firstrecstep) firstrecstep=False;
      }  /* End of loop over time steps */
  
      /*printf("Error: x: %lg ,y: %lg ,z: %lg\n",recoil.ddx,recoil.ddy,recoil.ddz);*/

       /*spread the deposited energy, moves values from temp-table to depen_table*/
      /*if no spreading, the values are already in depen_table*/
      if(physical.Damage)
	if(physical.spreadFD)
	  {
	  r2=recoil.sz;  
	  if (reccalc.Trange == 1) r2=recoil.r_proj;
	  if (reccalc.poly) r2+=poly.range;	 
          FD=(recoil.eprev-recoil.energy)-elstop.Se*recoil.dr-elstop.totalfirsov*recoil.dr;
          spread_fdn(r2,box.movelim.z,&recresult,FD,3);       
	  }

      if(physical.sputtering)
	{
        check_depencalc(&physical,&elstop,&recoil,&reccalc,&box,&poly,&recresult,&layers,1,&Ymean);
	printf("Mean Y after %d ions: %g ,sputtered surface: %g \n",nrec,Ymean/nrec,physical.sput_surf);
	}

      /*if(t!=0) printf("Ion flux density=%g [ions/cm^2/s]\n",1*physical.Dosescaling/((5.43e-8)*(5.32e-8)*t));*/
      endrec(at,&recoil,&reccalc,&recresult,&poly,&file,&nat,nrec,recerror,recbound,&physical,issput);
      DEBUGSRR("Recoil ended with",recoil.r_proj,recoil.sz); DEBUGS("\n");
      if (!reccalc.poly || !poly.boundreached) {
	 if (debug) eltime();
	 noutput++;
	 if (noutput>=file.output) { 
	    noutput=0; 
	    printresults(&file,&reccalc,&recresult,&unit,&physical,&box,&elstop,nrec);
	    if (reccalc.Estat) getestat(&recoil,&reccalc,&file,True); 
	    oflush();
	 }
	 nrec++;
      }
      if (reccalc.collisionsout>=0.0) fflush(file.fpcollisions);

   } /* End of loop over recoil events */

   if(file.movie) fclose(fmov);

   if (reccalc.collisionsout>=0.0) close(file.fpcollisions);

   nrec--;
   fprintf(file.fpo,"\nEnd of simulation of %d events at ",nrec); eltime();
   printresults(&file,&reccalc,&recresult,&unit,&physical,&box,&elstop,nrec);
   if (reccalc.Estat) getestat(&recoil,&reccalc,&file,True); 

   
} /* END OF PROGRAM */


