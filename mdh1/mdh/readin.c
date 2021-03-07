
#include "common.h"

#define RIC(name,val,mess) readic(#name ":",&(name),(int) (val),mess,file->fppa,file->fps)
#define RRC(name,val,mess) readdc(#name ":",&(name),(real) (val),mess,file->fppa,file->fps)
#define RDC(name,val,mess) readdc(#name ":",&(name),(double) (val),mess,file->fppa,file->fps)
#define RSC(name,val,mess) readsc(#name ":",&(name),(char *) (val),mess,file->fppa,file->fps)

char *asplit(char name[], char ind[], int i);
char splits[120];
int readstat;

#define RICA(name,ind,val,mess) { strncpy(splits,asplit(#name ":",#ind,ind),120); \
				  readstat=readic(splits,(&(name)),(int) (val),mess, \
						  file->fppa,file->fps); \
			       }
#define RRCA(name,ind,val,mess) { strncpy(splits,asplit(#name ":",#ind,ind),120); \
				  readstat=readdc(splits,(&(name)),(real) (val),mess, \
						  file->fppa,file->fps); \
			       }
#define RDCA(name,ind,val,mess) { strncpy(splits,asplit(#name ":",#ind,ind),120); \
				  readstat=readdc(splits,(&(name)),(double) (val),mess, \
						  file->fppa,file->fps); \
			       }
#define RSCA(name,ind,val,mess) { strncpy(splits,asplit(#name ":",#ind,ind),120); \
				  readstat=readsc(splits,(&(name)),(char *) (val),mess, \
						  file->fppa,file->fps); \
			       }

#define RICM(name,ind1,ind2,val,mess) { strncpy(splits,asplit(#name ":",#ind1,ind1),120); \
				        strncpy(splits,asplit(splits,#ind2,ind2),120); \
				        readstat=readic(splits,(&(name)),(int) (val),mess, \
						  file->fppa,file->fps); \
			       }

#define RDCM(name,ind1,ind2,val,mess) { strncpy(splits,asplit(#name ":",#ind1,ind1),120); \
				        strncpy(splits,asplit(splits,#ind2,ind2),120); \
				        readstat=readdc(splits,(&(name)),(double) (val),mess, \
						  file->fppa,file->fps); \
			       }

#define RDCM3(name,ind1,ind2,ind3,val,mess) { strncpy(splits,asplit(#name ":",#ind1,ind1),120); \
				              strncpy(splits,asplit(splits,#ind2,ind2),120); \
					      strncpy(splits,asplit(splits,#ind3,ind3),120); \
  				              readstat=readdc(splits,(&(name)),(double) (val),mess,\
						  file->fppa,file->fps); \
			       }

readin(type,time,layers,potcrit,pot,recoil,reccalc,box,physical,file,gen,elstop)
struct type *type;
struct time *time;
struct layers *layers;
struct potcrit *potcrit;
struct pot *pot;
struct recoil *recoil;
struct reccalc *reccalc;
struct box *box;
struct physical *physical;
struct file *file;
struct gen *gen;
struct elstop *elstop;
{
   static FILE *a,*b;
   int i=0;
   char s[120];
   
   a=fopen(file->param,"r");
   b=fopen(file->startf,"w");
   
   file->fppa=a;
   file->fps=b;

   if (maxdebug) puts("Going into time read");
   readtime(file,time); oflush();
   if (maxdebug) puts("Going into type read");
   readtype(file,type,reccalc); oflush();
   if (maxdebug) puts("Going into layers read");
   readlayers(file,layers); oflush();
   if (maxdebug) puts("Going into potential read");
   readpot(file,pot,type);
   if (maxdebug) puts("Going into box size read");
   readbox(file,box);
   if (maxdebug) puts("Going into recoil calculation read");
   readreccalc(file,reccalc,box);
   if (maxdebug) puts("Going into potcrit read");
   readpotcrit(file,potcrit,box);
   if (maxdebug) puts("Going into physical data read");
   readphysical(file,physical);
   if (maxdebug) puts("Going into miscellaneous data read");
   readgen(file,gen);
   if (maxdebug) puts("Going into file read");
   readfile(file);
   if (maxdebug) puts("Going into elstop param read");
   readelstopparam(file,elstop,reccalc);
   
   fflush(file->fps); oflush();

}
/*--------------------------------------------------------------------------*/
readtime(file,time)
struct file *file;
struct time *time;
{

   fprintf(file->fps,"\nTIME PARAMETERS\n\n");
   RRC(time->Dt,0.1,"Initial time step (fs) (NOT USED)");
   RRC(time->Dtmax,2.0,"Maximum time step (fs)");
   RRC(time->Kt,0.1,"Time step increasing criterion (Å)");
   RRC(time->Et,300.0,"Time step increasing criterion in collisions (eV)");
   RRC(time->Dtini,2.0,"Time step in inistate calc. (fs)");
   RRC(time->Timeini,200.0,"Length of inistate calc. (fs)");

}
/*--------------------------------------------------------------------------*/
readtype(file,type,reccalc)
struct file *file;
struct type *type;
struct reccalc *reccalc;
{
   int i;
   real psum=0.0;

   fprintf(file->fps,"\nTYPE PARAMETERS\n\n");
   reccalc->change=False;
   oflush(stdout);
   for (i=0;i<MAXTYPE;i++) {
      RICA(type[i].Z,i,14,     "Atomic number ");
      RRCA(type[i].m,i,28.0855,"Atom mass     ");
      RRCA(type[i].q,i,0.0,    "Atom charge   ");
      RRCA(type[i].prob,i,1.0, "Atom probability of occuring, if 1 no effect");
      if (type[i].prob != 1.0 && type[i].prob != 0.0) {
	 psum+=type[i].prob; 
         reccalc->change=True; 
      }
      fprintf(file->fps,"\n");
   }
   if (reccalc->change) {
      fprintf(file->fpo,"Atom change calculations performed\n"); oflush();
      if (oabs(psum-1.0)>1e-4) {
	 fprintf(file->fpo,"Warning: probability sums not = 1.0: psum=%g;\n",psum);
	 fprintf(file->fpo," are you sure you know what you're doing ??\n");
      }
   }
}
/*--------------------------------------------------------------------------*/
readlayers(file,layers)
/*
  See subroutine getlayer in borders.c
*/
struct file *file;
struct layers *layers;
{
   int i;
   static struct layer *tmp;
   
   fprintf(file->fps,"\nLAYERS PARAMETERS\n\n"); oflush();
   RIC(layers->Nlayers,1,"Number of different layers");
   MAXLAYER=layers->Nlayers;
   tmp=(struct layer *) malloc(MAXLAYER*sizeof(struct layer));
   layers->layer=tmp;

   DEBUG(MAXLAYER); oflush();
   if (layers->Nlayers==1) layers->onelayer=True;
   else layers->onelayer=False;

   RIC(layers->Ninistate,1,"Number of inistate calcs for each layer");
   MAXINISTATE=layers->Ninistate;
   DEBUG(MAXINISTATE);

   RIC(layers->Natomsmax,300,"Maximum number of atoms in one layer");
   MAXATOM=layers->Natomsmax;
   DEBUG(MAXATOM);

   for (i=0; i< layers->Nlayers && i<MAXLAYER;i++) {
      fprintf(file->fps,"Layer number %d:\n",i);
      RRCA(layers->layer[i].minz,i,-1e9,   "Minimum Z  ");
      RRCA(layers->layer[i].maxz,i, 1e9,   "Maximum Z  ");
      RRCA(layers->layer[i].density,i,2.33,"Density (not used)");
      RDCA(layers->layer[i].debyet,i,0, "Debye temp (0 if N/A)");
      RICA(layers->layer[i].ltype,i,0,     "Layer type, i.e. index > 0 of coord file name");
      RRCA(layers->layer[i].Timeinior,i,-1.0, "Override for time of inistate calc.");
      RRCA(layers->layer[i].Probability,i,1.0, "Probability for layer selection");
      RICA(layers->layer[i].Probmode,i,0, "Mode of probability for layer selection");
      RRCA(layers->layer[i].Thickness,i,0.0, "Thickness of selected layer for Probmode=2");
      if (layers->layer[i].Probability>1.0) 
	 fprintf(fpdebug,"readlayers: Warning !! : layer probability > 1\n");
      fprintf(file->fps,"\n");
   }
   oflush();
}

/*--------------------------------------------------------------------------*/
readpot(file,pot,type)
struct file *file;
struct pot *pot;
struct type *type;
{
   int i,j,k;
   double a0=0.529177;    /* Bohr atom radius */
   logical indflag=False;
   static char *attrtype[] = {
         "Zero: unknown", "Mazzone Silicon potential","Morse"
      };

/* REPULSIVE POTENTIAL */

   fprintf(file->fps,"\nREPULSIVE POTENTIAL PARAMETERS\n\n");

   fprintf(file->fps,"\nPotential type for each pair: \n");
   fprintf(file->fps,"1 = ZBL type-potential using read in params \n");
   fprintf(file->fps,"2 = Potential and derivative read in from file and spline interpolated\n");

   for (k=0;k<MAXTYPE;k++) {
      for (j=0;j<MAXTYPE;j++) {
	 RICM(pot->rep.type[k][j],k,j,1," ");
    	 if (pot->rep.type[k][j] == 1) ;
	 else if (pot->rep.type[k][j] == 2) ;
	 else if (pot->rep.type[k][j] == 3) ;
	 else {
	    fprintf(file->fpdb,"readpot ERROR: Unknown repulsive potential type %d fopr pair %d %d\n",
		    pot->rep.type[k][j],k,j);
	    exit(0);
	 }
      }
   }

/* First find the default parameters for all interactions */

   for (i=0;i<MAXPOTEXP;i++) {
      RDCA(pot->rep.ai[i],i,0.0,"   ");
      if (pot->rep.ai[i]==0.0) break;
      RDCA(pot->rep.bi[i],i,0.0,"   ");
   }
   pot->rep.nterms=i;
   for (;i<MAXPOTEXP;i++) { pot->rep.ai[i]=0.0; pot->rep.bi[i]=0.0; }

   if (pot->rep.ai[0]==0.0) {
      fprintf(file->fps,"ai[0] is 0.0 => using ZBL values\n");
      pot->rep.ai[0]=0.1818;  pot->rep.bi[0]=-3.2;
      pot->rep.ai[1]=0.5099;  pot->rep.bi[1]=-0.9423;
      pot->rep.ai[2]=0.2802;  pot->rep.bi[2]=-0.4029;
      pot->rep.ai[3]=0.02817; pot->rep.bi[3]=-0.2016;
      pot->rep.nterms=4; 
   }
   fprintf(file->fps,"\nDefault repulsive potential factors (used for all interactions !): \n\n");
   for(i=0;i<MAXPOTEXP;i++) {
      fprintf(file->fps,"%d %.10lg %.10lg\n",i,pot->rep.ai[i],pot->rep.bi[i]);
   }

   MAXDEBUG((real) pot->rep.nterms);

/* Then look for exponential params for specific interactions */

   fprintf(file->fps,"\nIndividual repulsive potential factors (used for specific atom pair):\n\n");
   for (k=0;k<MAXTYPE;k++) {
      for (j=0;j<MAXTYPE;j++) {
	 for (i=0;i<MAXPOTEXP;i++) {
	    RDCM3(pot->rep.a[k][j][i],k,j,i,0.0,"\0");
	    if (pot->rep.a[k][j][i]==0.0) break;
	    RDCM3(pot->rep.b[k][j][i],k,j,i,0.0,"\0");
	    pot->rep.a[j][k][i]=pot->rep.a[k][j][i];
	    pot->rep.b[j][k][i]=pot->rep.b[k][j][i];
	    if (i==0) {
	       indflag=True;
	       fprintf(file->fpo,"Found individual repulsive potential between atom types %d and %d\n",k,j);
	    }
	    fprintf(file->fps,"Type1 %d Type2 %d i %d: a %.10lg b %.10lg\n",k,j,i,pot->rep.a[k][j][i],pot->rep.b[k][j][i]);
	 }
      }
   }
   if (!indflag) fprintf(file->fps,"None found !\n");

/* Finally overwrite those atom pairs that don't have an individual potential  */ 

   pot->rep.k_e= 14.39976;
   fprintf(file->fps,"pot->rep.k_e %g\n",pot->rep.k_e);
   for (i=0;i<MAXTYPE;i++) {
      for (j=0;j<MAXTYPE;j++) {
	 if (pot->rep.a[i][j][0] == 0.0) {
	    for (k=0;k<MAXPOTEXP;k++) {
	       pot->rep.a[i][j][k]=pot->rep.ai[k];
	       pot->rep.b[i][j][k]=pot->rep.bi[k];
	       MAXDEBUGSRRRRR("  ",i,j,k,pot->rep.a[i][j][k],pot->rep.b[i][j][k]);
	    }
	 }
         pot->rep.fact[i][j]=pot->rep.k_e*type[i].Z*type[j].Z;
         pot->rep.au[i][j]=0.8854*a0/(pow((double)type[i].Z,0.23)+pow((double)type[j].Z,0.23));
	 fprintf(file->fps,"typei %d typej %d fact %g au %g\n",i,j,pot->rep.fact[i][j],pot->rep.au[i][j]);
      }
   } 





   RRC(pot->rep.scale,1.0,"\nRepulsive potential scaling factor");
   
/* Attractive potential */

   fprintf(file->fps,"\nATTRACTIVE POTENTIAL PARAMETERS\n\n");

   RIC(pot->attr.type,1,"Attractive potential type");
   fprintf(file->fps,"Type means %s\n",attrtype[pot->attr.type]); 

   if (pot->attr.type == 1) {
      RRC(pot->attr.mazD,4.5,"   D (eV)");
      RRC(pot->attr.mazalpha,0.2,"   alpha (1/Å)");
      RRC(pot->attr.mazr1,2.35,"   Nearest neighbour dist (Å)");
      RRC(pot->attr.mazK,1.8,"   Harmonic const (eV/Å)");
      RRC(pot->attr.mazr2,3.84,"   Second nearest neighbourdist (Å)");
      RRC(pot->attr.mazd,0.33,"    Width of harmonic well (Å)");
   }

   if (pot->attr.type == 2) {
      RIC(pot->attr.morsetype,0,"Morse potential type (if 0 values readin)");
      if (pot->attr.morsetype==0) {
         RRC(pot->attr.morsealpha[0],1.0,"Morse alpha");
         RRC(pot->attr.morser1[0],3.0,"Morse r1");
         RRC(pot->attr.morseD[0],0.5,"Morse D");
      }
   } 


}
/*--------------------------------------------------------------------------*/
readbox(file,box)
struct file *file;
struct box *box;
{
   fprintf(file->fps,"\nSIMULATION CELL PARAMETERS\n\n");

   RRC(box->Min.x,0.0,"Box limit ");
   RRC(box->Min.y,0.0,"Box limit ");
   RRC(box->Min.z,0.0,"Box limit ");
   RRC(box->Max.x,10.86,"Box limit ");
   RRC(box->Max.y,10.86,"Box limit ");
   RRC(box->Max.z,10.86,"Box limit ");
   box->size.x=box->Max.x-box->Min.x;
   box->size.y=box->Max.y-box->Min.y;
   box->size.z=box->Max.z-box->Min.z;
   RRC(box->movelim.x,box->size.x/4.0,"Limit for moving cell");
   RRC(box->movelim.y,box->size.y/4.0,"Limit for moving cell");
   RRC(box->movelim.z,box->size.z/4.0,"Limit for moving cell");
   RIC(box->Iniscaleint,50,"Time step interval between cm scalings in inistate calc.");
}
/*--------------------------------------------------------------------------*/
readreccalc(file,reccalc,box)
struct file *file;
struct reccalc *reccalc;
struct box *box;
{
   fprintf(file->fps,"\nRECOIL CALCULATION PARAMETERS\n\n");
   RIC(reccalc->type,1,"Recoil calculation type");
   RIC(reccalc->Atype,1,"Recoil particle type");
   RIC(reccalc->Ncalc,10000,"Number of recoils to be calculated");
   RIC(reccalc->Trange,1,"Type of range used for result (0 sz 1 proj 2 chord 3 path)");

   RRC(reccalc->E0,10000.0,"Original energy in eV");
   RRC(reccalc->Estatzmax,5000.0,"Maximum z for arrays in Estat calc");
   RRC(reccalc->Emin,1.0,"Minimum energy before recoil is stopped");
   RRC(reccalc->Zmax,1e30,"Z at which to end recoil calculation");

   RRC(reccalc->Startmin.x,box->Min.x+box->movelim.x,"Starting position min limit ");
   RRC(reccalc->Startmin.y,box->Min.y+box->movelim.y,"Starting position min limit ");
   RRC(reccalc->Startmin.z,-box->movelim.z,"Starting position min limit ");
   RRC(reccalc->Startmax.x,box->Max.x-box->movelim.x,"Starting position max limit ");
   RRC(reccalc->Startmax.y,box->Max.y-box->movelim.y,"Starting position max limit ");
   RRC(reccalc->Startmax.z,-box->movelim.z,"Starting position max limit ");

   RRC(reccalc->Theta0,0.0,"Velocity direction (degrees)");
   RRC(reccalc->Fii0,0.0,"Velocity direction (degrees)");
   RRC(reccalc->Thetamax,0.0,"Velocity direction limit (degrees)");
   RRC(reccalc->Fiimax,0.0,"Velocity direction limit (degrees)");
   reccalc->Theta0*=3.1415927/180.0;
   reccalc->Thetamax*=3.1415927/180.0;
   reccalc->Fii0*=3.1415927/180.0;
   reccalc->Fiimax*=3.1415927/180.0;

   RRC(reccalc->Polysize,0.0,"Polycrystalline domain size");
   RRC(reccalc->Polysizesigma,reccalc->Polysize/5.0,"Poly domain size sigma");
   RRC(reccalc->Polydeltaangle,0.0,"Maximum change in direction, in degrees");
   reccalc->Polydeltaangle*=3.1415927/180.0;
   if (reccalc->Polysize==0.0) reccalc->poly=False;
   else reccalc->poly=True;

   RIC(reccalc->Forceseed,0,"FORCED SEED FOR RECOIL CALC, DEBUG ONLY !!");
   RIC(reccalc->Outputstep,0,"Interval for coordinate output");

   RRC(reccalc->recspeczmin,-1e30,"Minimum z for primary recoil spectrum calc.");
   RRC(reccalc->recspeczmax,1e30,"Maximum z for primary recoil spectrum calc.");

}
/*--------------------------------------------------------------------------*/
readpotcrit(file,potcrit,box)
struct file *file;
struct potcrit *potcrit;
struct box *box;
{
   fprintf(file->fps,"\nPOTCRIT PARAMETERS\n\n");
   RRC(potcrit->R0rec,3.0,"Potential cutoff radius in recoil calc.        ");
   RRC(potcrit->Rmrec,4.0,"Neighhbour list cutoff radius in recoil calc.  ");
   RRC(potcrit->R0ini,4.3,"Potential cutoff radius in inistate calc.      ");
   RRC(potcrit->Rmini,box->size.x/2.0-0.1,"Neighhbour list cutoff radius in inistate calc.");


}
/*--------------------------------------------------------------------------*/
readphysical(file,physical)
struct file *file;
struct physical *physical;
{
   fprintf(file->fps,"\nPHYSICAL SIMULATION PARAMETERS\n\n");

   RRC(physical->Tini,300.01,"Required temperature in inistate calc.");
   RRC(physical->Trec,physical->Tini,"Required temperature in recoil calc.");
   
}
/*--------------------------------------------------------------------------*/
readgen(file,gen)
struct file *file;
struct gen *gen;
{
   fprintf(file->fps,"\nMISCELLANEOUS GENERAL PARAMETERS\n\n");

   RIC(gen->seed,438724623,"Seed number for random number generator");

}
/*--------------------------------------------------------------------------*/
readfile(file)
struct file *file;
{
   fprintf(file->fps,"\nFILE HANDLING PARAMETERS\n\n");

   RIC(file->output,40,"Interval between result output in range calc.");

}
/*--------------------------------------------------------------------------*/
readelstopparam(file,elstop,reccalc)
struct file *file;
struct elstop *elstop;
struct reccalc *reccalc;
{
   fprintf(file->fps,"\nELSTOP PARAMETERS\n\n");
   RIC(elstop->mode,0,"El. stopping mode: if 1 include straggling");

   RRC(elstop->scale,1.0,"Elstop modifying factor, see elstop.c");
   RRC(elstop->dscale,0.0,"Linear elstop modifying factor");
   RRC(elstop->sqscale,0.0,"Quadratic elstop modifying factor");

   if (elstop->dscale!=0.0) {
      fprintf(file->fpo,"\nAdvanced elstop scaling selected:\n");
      fprintf(file->fpo,"Elstop at zero energy will be scaled %g\n",elstop->scale);
      fprintf(file->fpo,"at recoil energy %g scaling is %g\n\n",reccalc->E0,
         elstop->scale+reccalc->E0*(elstop->dscale+elstop->sqscale*reccalc->E0));
   }

}
/*--------------------------------------------------------------------------*/


