#define MAIN 

#include "common.h" 

#define VERSION "mdh V3. REED Aug 30 2002"

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
  3.0Zeta  rewritten moveslice, added debuggin option -pr z, usefull for monitoring the surrounding
           of the ion at depth z.
  3.0REED  added reed-option, -dens option for calculating el. dens. and Se spectra, 
           NOTE! use this model only for calculating range-profiles! Other data calculation is not tested, it may work or not..
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
   logical recerror;

   int i;       /* way of stopping */

  
   real Ymean=0.0,ddz=5.43,r2,FD,expath=0;
   int j,k;
   logical templog;
   struct atom *attemp;
   struct atomflag *atflagtemp;
   struct atomex *atextemp;
   point origmovelim,origboxsize;
   point tat;
   FILE *fad,*fmov;
   int frame=0;
  

   /*reed */
   int nr,nreed=0,nats[MAXSPLITATOMS];
   struct recoil reedrecoils[MAXSPLITATOMS];
   struct box reedbm[MAXSPLITATOMS];
   static struct atom *splitat;          /* Atom coordinates */
   static struct atomex *splitatex;      /* Miscellaneous atom properties */
   static struct atomflag *splitatflag;  /* flags etc. */


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

   
   if(physical.Damage || physical.sputtering)
     {
     physical.Dosescaling=physical.Dose*(1e-8)*(1e-8)*(5.43)*(5.43)/reccalc.Ncalc;
     printf("DOSE SCALING in damage buildup model=%g\n",physical.Dosescaling);
     }
  

/*large atom tables for the simulation cells*/
   at= (struct atom *) malloc(MAXATOM*sizeof(struct atom));
   attemp = (struct atom *) malloc(MAXATOM*sizeof(struct atom)); 
   atex= (struct atomex *) malloc(MAXATOM*sizeof(struct atomex));
   atextemp= (struct atomex *) malloc(MAXATOM*sizeof(struct atomex)); 
   atflag= (struct atomflag *) malloc(MAXATOM*sizeof(struct atomflag));
   atflagtemp= (struct atomflag *) malloc(MAXATOM*sizeof(struct atomflag)); 
   /*reed*/
   if(reccalc.reed)
     {
     splitatex= (struct atomex *) malloc(MAXSPLITATOMS*MAXATOM*sizeof(struct atomex)); 
     splitatflag= (struct atomflag *) malloc(MAXSPLITATOMS*MAXATOM*sizeof(struct atomflag));
     splitat= (struct atom *) malloc(MAXSPLITATOMS*MAXATOM*sizeof(struct atom));
     }

   iniat= (struct iniat *) malloc(MAXLAYER*MAXAMORPH*MAXINISTATE*MAXATOM*sizeof(struct iniat));
   MAXDEBUGSR("Initialized iniat array to ",MAXLAYER*MAXAMORPH*MAXINISTATE*MAXATOM*sizeof(struct iniat));

   if (debug) puts("Going into parameter initialization");
   initparams(type,&unit,&time,&gen,&file);

   if (debug) puts("Going into atom coordinate reading");
   readcoords(iniat,&layers,&file,type,&box,&unit,&physical);



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
   
   origmovelim=box.movelim;  
   origboxsize=box.size;

   
   if(file.movie) fmov=fopen("at.movie","w");

   if(reccalc.reed) for(i=0;i<MAXSPLITATOMS;i++) recresult.splitborders[i]=1000000; /*huge default values */

   for (nrec=1;nrec<=reccalc.Ncalc;) {
     /*reed*/
     if(reccalc.reed)
       for(j=0; j<MAXSPLITATOMS; j++)
          zeroatom(splitat+j,splitatex+j,splitatflag+j);	

     nr=1; /*virtual atom counter for each recoil ion*/
     nreed=0; /*recoil atom counter, goes up to the total number of simulated ions for each real ion, 0 is the original ion, 1->nr are the divided ones.*/
     reedrecoils[nreed]=recoil;

    do{  /*reed loop	      */
      firstrecstep=True;
      if(nreed==0)
	{
        getlayer(&layers,reccalc.Startmax.z,&nlayer,firstrecstep);
	getamorphstate(&namorph,&recresult,reccalc.Startmax.z,box.movelim.z,&physical,&layers,nlayer); 
	}
      else
	{	  
	  getlayer(&layers,reedrecoils[nreed].sz,&nlayer,firstrecstep);
	  getamorphstate(&namorph,&recresult,reedrecoils[nreed].sz,box.movelim.z,&physical,&layers,nlayer); 
	}
      
      if(physical.scaling)
	{
	  /*if scaling done in the readcoords.c, the movelims must be set right for the next ion.*/
	if(debug) fprintf(fpdebug,"Bscaling is set on and movelim.sizes are scaled for the first layer\n");
	box.movelim=origmovelim;box.size=origboxsize;
	}
      
      
      if(nreed==0) nat=layers.layer[nlayer].natoms[namorph];
      else nat=nats[nreed];
      DEBUGSRRR("Starting recoil ",nrec,nat,poly.ndomain);
      /* Zero arrays, set open boundaries ... */ 
      if(nreed==0)
	{
         startrec(iniat,at,atex,atflag,&layers,&nlayer,&ninistate,&box,&potcrit,&time,
	       &reccalc,&recresult,&poly,&elstop,&namorph);  
	 /*printf("BOX(%g,%g,%g)\n",box.movement.x,box.movement.y,box.movement.z);*/
	}
      else{
	  ninistate=unirand(0)*layers.Ninistate;
	  if(ninistate>=layers.Ninistate) ninistate=layers.Ninistate-1;
	  /*copyinilayer(iniat,at,atex,atflag,nlayer,ninistate,namorph);*/
	  for(i=0; i<nat; i++)
	     copyatom(&splitat[nreed*MAXATOM+i],&splitatex[nreed*MAXATOM+i],&splitatflag[nreed*MAXATOM+i],
		      &at[i],&atex[i],&atflag[i]);
	  box=reedbm[nreed];
	  /*printf("BOX(%g,%g,%g)\n",box.movement.x,box.movement.y,box.movement.z);	  */
	  }


      /*check if first layer has different dimensions than simulation box=at[]*/
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

      if(nreed==0)
	{
        createrecatom(&reccalc,&recoil,&maxima,&nat,at,atflag,&poly,nrec);
	/*printf("\nREED: new recoil atom nr.%d (z;%g),weight %g\n",nrec,reccalc.Startmax.z,pow(0.5,recoil.splitnum));*/
	}
      else 
	{
	  /*copyatom(&splitat[nreed],&splitatex[nreed],&splitatflag[nreed],&at[nat],&atex[nat],&atflag[nat]);*/
	  /*nat++;*/
	recoil=reedrecoils[nreed];
	getmaximumprop(&maxima,at,nat); 
	/*printf("\n   REED: old recoil(%d), virtual atom nr.%d / %d (z:%g),weight %g\n",nrec,nreed,nr-1,reedrecoils[nreed].sz,pow(0.5,recoil.splitnum));*/
	/*printf("recoilnumber splitted:%d\n",recoil.number);*/
	/*printf("recoil energy for split :%g\n",recoil.energy);*/
	}
      
      if(physical.sput_surf>0 && nreed==0)  /* surface modifications*/
        erodesurface(at,atex,atflag,attemp,atextemp,atflagtemp,&recoil,&physical,&box,&nat);

      if(nreed==0)
      for(j=0; j<ZBOXES; j++)
       recresult.depennucl_boxtemp[j]=0.0;

      /* Select neigbour list mode, see accelerations.c */
      if (!file.fullcalc) ngbrlistmode=1;   
      else ngbrlistmode=0; 

      time.dtreal=time.Dt*1e-15;
      inittsteps(&time,&unit);
      calctsteps(&time,&unit,&maxima,&recoil,True);
      
      if(nreed==0)
      for(i=0;i<layers.Nlayers;i++) /*setting box movement to zero for all layers.*/
	for(k=0;k<layers.layer[i].Namorph; k++)
	  {
	  layers.layer[i].front[k].x=0.0;
          layers.layer[i].front[k].y=0.0;
          layers.layer[i].front[k].z=0.0;
	  }

      oflush();
      update=True; recerror=False;
      while( (recoil.energy>reccalc.Emin || recoil.incollision) 
	    && !(reccalc.type!=2 && recoil.range<(2.0*reccalc.Startmin.z)) ) {

	/*reed*/
	if(reccalc.reed)
	if(splitlimitcheck(&recoil,&recresult))
	  {
	  if(nr<MAXSPLITATOMS) 
	    {
	    for(i=0; i<nat; i++)
	       copyatom(&at[i],&atex[i],&atflag[i],
			&splitat[nr*MAXATOM+i],&splitatex[nr*MAXATOM+i],&splitatflag[nr*MAXATOM+i]);
	    recoil.splitnum++;
	    reedrecoils[nr]=recoil;
	    nats[nr]=nat;
	    reedbm[nr]=box;
	  
	    /* if(nreed==0)
	      printf("REED: splitted real atom nr.%d (z:%g) to virtual atom nr %d,weight %g\n",nrec,recoil.sz,nr,pow(0.5,recoil.splitnum));
	    else 
	      printf("REED: splitted virtual atom nr.%d (z:%g) to virtual atom nr %d,weight %g\n",nreed,recoil.sz,nr,pow(0.5,recoil.splitnum));
	    printf("recoilnumber splitting:%d\n",recoil.number);
	    printf("recoil energy when splitting :%g\n",recoil.energy);*/
	    nr++;
	    }
	  }
	
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
		   file.recoilspec,&physical);  /*JARKKO-damage*/
	
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
	 if(recoil.sz>=physical.sput_surf)  /*JARKKO*/
           subelstop(at,&time,&recoil,&elstop,&unit,&type,layers.layer[nlayer],&potcrit,nlayer,&file); 
	 else{elstop.totalfirsov=0;
	      elstop.Se=0;}

  
	 if (reccalc.Estat) getestat(&recoil,&reccalc,&file,False);


	 check_depencalc(&physical,&elstop,&recoil,&reccalc,&box,&poly,&recresult,&layers,0,&Ymean);
	
	 if (file.printdepen || file.printstop || file.printavvel) 
	    calcmisc(&recoil,&elstop,&reccalc,&recresult,&poly,&file,
		      firstrecstep,t/unit.time[recoil.attype],&physical,nreed);  
	                   /* ,layers.layer[nlayer].ndensity); */
	 MAXDEBUGS("Getting new t step\n");
         calctsteps(&time,&unit,&maxima,&recoil,False);
	 
	 if (reccalc.poly && recoil.r_proj > poly.cursize) {
	    exitpolydomain(&poly,&recoil);
	    break;   /* Exit loop over time steps */
	 }

	 if (recoil.energy > 1.3* reccalc.E0) {
	    fprintf(fpdebug,"Warning: Recoil n. %d has insane energy %g. Discarding it... try smaller timestep!\n",nrec,recoil.energy);
	    printrec(&recoil,t,0);
	    printcoint(at,atflag,nat,&box,t,0);
	    recerror=True;   /* Discard this recoil event */
	    break;   /* Exit loop over time steps */
	 }

	 if (recoil.range>reccalc.Zmax || recoil.range>layers.layer[layers.Nlayers-1].maxz) break;
	 if (recoil.sz<(2.0*reccalc.Startmin.z)) break; 

         if (file.printcoords) printcoint(at,atflag,nat,&box,t,reccalc.Outputstep);
         if (file.printrec || maxdebug) printrec(&recoil,t,reccalc.Outputstep);
         t+=time.dtreal;
	 if (firstrecstep) firstrecstep=False;
      }  
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
	  if(reccalc.reed) FD*=pow(0.5,recoil.splitnum);
          spread_fdn(r2,box.movelim.z,&recresult,FD,3);       
	  }

      if(physical.sputtering)
	{
        check_depencalc(&physical,&elstop,&recoil,&reccalc,&box,&poly,&recresult,&layers,1,&Ymean);
	printf("Mean Y after %d ions: %g ,sputtered surface: %g \n",nrec,Ymean/nrec,physical.sput_surf);
	}

      /*if(t!=0) printf("Ion flux density=%g [ions/cm^2/s]\n",1*physical.Dosescaling/((5.43e-8)*(5.32e-8)*t));*/
      endrec(at,&recoil,&reccalc,&recresult,&poly,&file,&nat,nrec,recerror,&physical,nreed);
      DEBUGSRR("Recoil ended with",recoil.r_proj,recoil.sz); DEBUGS("\n");

      /*printf("End coords stopped: %g,%g,%g ,%g,%g,%g\n",recoil.sx,recoil.sy,recoil.sz,recoil.vx,recoil.vy,recoil.vz);*/

       nreed++;
      }while(nreed<nr); /*end of reed ions*/
      if(reccalc.reed) splitlimitcalc(&recresult); /*new borders*/

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

   } /* End of loop over recoil events */

   if(file.movie) fclose(fmov);

   nrec--;
   fprintf(file.fpo,"\nEnd of simulation of %d events at ",nrec); eltime();
   printresults(&file,&reccalc,&recresult,&unit,&physical,&box,&elstop,nrec);
   if (reccalc.Estat) getestat(&recoil,&reccalc,&file,True); 

   
} /* END OF PROGRAM */










