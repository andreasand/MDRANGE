#include "common.h"

inistatecalc(iniat,pot,potcrit,time,unit,type,file,box,physical,layers)
struct iniat *iniat;
struct pot *pot;
struct potcrit *potcrit;
struct time *time;
struct unit *unit;
struct type *type;
struct file *file;
struct box *box;
struct physical *physical;
struct layers *layers;
{
   struct atom *at;         
   struct atomex *atex;
   struct atomflag *atflag;
   
   struct box box2;
   
   int nlayer,ninistate,nat,namorph;
   real t,T,tmax;
   int i,intcounter;
   real deltarmean;
   point *origpos;
   logical update=False;
   double sd;
   int j,k,n;

   at = (struct atom *) malloc(MAXATOM*sizeof(struct atom));
   atex= (struct atomex *) malloc(MAXATOM*sizeof(struct atomex));
   atflag= (struct atomflag *) malloc(MAXATOM*sizeof(struct atomflag));
   origpos= (point *) malloc(MAXATOM*sizeof(point));

   gaussianvels(iniat,type,unit,layers,physical->Tini);
   fprintf(file->fpo,"\nStarting inistate calculations\n");
   fprintf(file->fpo,  "------------------------------\n\n");
   
   for (nlayer=0;nlayer<layers->Nlayers;nlayer++)
     {
    for(namorph=0; namorph<layers->layer[nlayer].Namorph; namorph++) {
     
      for (ninistate=0;ninistate<layers->Ninistate;ninistate++) {
         fprintf(file->fpo,"Layer # %d amorphstate # %d Inistate calc. # %d\n",nlayer,namorph,ninistate);

/*
   Do necessary initializations
*/

         time->dtreal=time->Dtini*1e-15;
         inittsteps(time,unit);
	 t=0.0;
	 nat=layers->layer[nlayer].natoms[namorph];

         copyinilayer(iniat,at,atex,atflag,nlayer,ninistate,namorph,nat);
         
         if (file->printcoordsi) printcoveac(at,nat,t);

	 DEBUG(T);
         box->periodic=True; box->warnoutside=True;
	 potcrit->r0=potcrit->R0ini;
	 potcrit->rm=potcrit->Rmini;
	 DEBUG(nat); 
	 for (i=0;i<nat;i++) { origpos[i].x=at[i].x; origpos[i].y=at[i].y; origpos[i].z=at[i].z; }
	 for (i=0;i<nat;i++) change(at,type,i);

/*
   Then the actual inistate calculation
*/

	 box2.periodic=box->periodic;  /*setting sizes for checking*/
	 box2.warnoutside=box->warnoutside;
	 box2.size=layers->layer[nlayer].boxsize[namorph];
 
	 intcounter=0;
	 tmax=time->Timeini;
	 if (layers->layer[nlayer].Timeinior >= 0.0) tmax=layers->layer[nlayer].Timeinior; 
	 DEBUGSRRR("Going into inistate loop",nlayer,ninistate,tmax);

         if(layers->layer[nlayer].debyet<=0) {  
	   printf("Debye T not set, attempting inistate MD simulation\n");

	    setcelltemp(at,atflag,nat,unit,2*physical->Tini); /* The factor 2 is not a bug */
	    getcelltemp(at,atflag,nat,unit,&T,&(physical->ekin));
	    
	    if (physical->Tini>0.0) {
	       for (t=0;t<=tmax;t+=time->Dtini) {
		  MAXDEBUGSR("\nStart of inistate loop\n",t);
		  solve1(at,time,atflag,box,nat);
		  borders(at,box2,nat);
		  accelerations(at,atex,atflag,pot,potcrit,box,physical,unit,nat,0,&update,&file);
		  solve2(at,time,atflag,box,nat);
		  getcelltemp(at,atflag,nat,unit,&T,&(physical->ekin));
		  MAXDEBUGSRRRR("ini",t,T,physical->ekin,physical->epot);
		  
		  intcounter++;
		  if (intcounter>=box->Iniscaleint) {
		     intcounter=0;
		     MAXDEBUGSR("Subtracting CM vel. and acceleration",t);
		     subtractcmvel(at,nat);
		     subtractcmacc(at,nat);
		  }
		  MAXDEBUG(T);
	       
		  if (file->printcoordsi) printcoveac(at,nat,t);
	       }
	    }
	    getmeandisplacement(&deltarmean,origpos,at,nat,box);
	    fprintf(file->fpo,"Mean displacement (sim.) %g\n\n",deltarmean);
	    setinilayer(iniat,at,nlayer,ninistate,namorph,nat); /*JARKKO-damage*/
       
	 }
	 else {
	    /* 21.3 is sqrt(9 hbar^2/k_B u) 
	       Expression comes from second order approximation of a nasty
	       elliptic integral
	    */
	    MAXDEBUGSRRR("Debye displ",physical->Tini,type[1].m,layers->layer[nlayer].debyet);

	   /* The debye/attr.pot switch parameter and the Debye temp have
           been added to struct layers in param.in */
	   for(i=0;i<nat;i++){
	     if (layers->layer[nlayer].debyet>0.0) {
	       sd=21.3*sqrt(physical->Tini/type[at[i].type].m)/(sqrt(3.0)*layers->layer[nlayer].debyet);
	     }
	     else {
	       fprintf(stderr,"Debye T not set or negative %g although needed %d\n",layers->layer[nlayer].debyet,nlayer);
	       exit(0);
	     }

	     at[i].x=at[i].x+sd*gaussianrand();
	     at[i].y=at[i].y+sd*gaussianrand();
	     at[i].z=at[i].z+sd*gaussianrand();
	     }
	 getmeandisplacement(&deltarmean,origpos,at,nat,box);
	 fprintf(file->fpo,"Mean displacement (Debye) %g\n\n",deltarmean);
	 setinilayer(iniat,at,nlayer,ninistate,namorph,nat);
	 }
       

      DEBUGSR("End of inistate calcs for layer ",nlayer);
       } /*inistate looppi*/
   
   DEBUGS("End of inistate calculations for all layers\n");
   } /* am-looppi*/
   DEBUGSR("End of inistate calcs for am state ",namorph); 
     }/*nlayer looppi*/
}

gaussianvels(iniat,type,unit,layers,T)
struct iniat *iniat;
struct type type[MAXTYPE];
struct unit *unit;
struct layers *layers;
real T;
{
  int nlayer,ninistate,n,namorph;
   real vx,vy,vz,vxsum=0.0,vysum=0.0,vzsum=0.0,v,vtot;
   real unirand(),gaussianrand();

   DEBUGS("Giving gaussian initial velocities to atoms\n");
  
   for (nlayer=0;nlayer<layers->Nlayers; nlayer++) {
     DEBUGSRR("Copying coordinates from layer to inistate arrays",nlayer,layers->Ninistate);
     for(namorph=0; namorph<layers->layer[nlayer].Namorph; namorph++) 
     {
       if(debug) printf("gaussianvel :layer %d, amstate %d\n",nlayer,namorph);
      for (ninistate=1;ninistate<layers->Ninistate; ninistate++) {
         for (n=0;n<layers->layer[nlayer].natoms[namorph]; n++) {
	   /*JARKKO-damage*/
            iniat[iniind(nlayer,namorph,ninistate,n)].x=iniat[iniind(nlayer,namorph,0,n)].x;
            iniat[iniind(nlayer,namorph,ninistate,n)].y=iniat[iniind(nlayer,namorph,0,n)].y;
            iniat[iniind(nlayer,namorph,ninistate,n)].z=iniat[iniind(nlayer,namorph,0,n)].z;
            iniat[iniind(nlayer,namorph,ninistate,n)].type=iniat[iniind(nlayer,namorph,0,n)].type;
	    /*printf("%g,%g,%g\n",iniat[iniind(nlayer,namorph,ninistate,n)].x,iniat[iniind(nlayer,namorph,ninistate,n)].y,iniat[iniind(nlayer,namorph,ninistate,n)].z);*/
         } /*atom*/
      } /*inistate*/
     } /*amorp*/
   } /* layer*/

   v=sqrt(2*T/unit->tempconv);
   if (debug) printf("Initial rms velocity %g, in m/s %g, in K %g\n",
		     v,v*unit->vel[0],T);
   for (nlayer=0;nlayer<layers->Nlayers;nlayer++) {
      MAXDEBUGSR("Number of layer ",nlayer);
      for(namorph=0; namorph<layers->layer[nlayer].Namorph; namorph++) 
      { 
      for (ninistate=0;ninistate<layers->Ninistate;ninistate++) {
	 MAXDEBUGSR("Number of inistate calc. in layer",ninistate);
	 vtot=0;
	 vxsum=vysum=vzsum=0.0; /*JARKKO-fix*/
         for (n=0;n<layers->layer[nlayer].natoms[namorph];n++) {
            vx=v*gaussianrand();
            vy=v*gaussianrand();
            vz=v*gaussianrand();
	    /*JARKKO-damage*/
            iniat[iniind(nlayer,namorph,ninistate,n)].vx=vx; vxsum+=vx;
            iniat[iniind(nlayer,namorph,ninistate,n)].vy=vy; vysum+=vy;
            iniat[iniind(nlayer,namorph,ninistate,n)].vz=vz; vzsum+=vz;
	    vtot+=sqrtf(vx*vx+vy*vy+vz*vz);
         }
	 vtot/=layers->layer[nlayer].natoms[namorph];
	 DEBUG(vtot);
         vxsum/=layers->layer[nlayer].natoms[namorph];
         vysum/=layers->layer[nlayer].natoms[namorph];
         vzsum/=layers->layer[nlayer].natoms[namorph];
	 if(debug) printf("vsums,initsate,amstate,layer : %g %g %g %d %d %d\n",vxsum,vysum,vzsum,ninistate,namorph,nlayer);
         for (n=0;n<layers->layer[nlayer].natoms[namorph];n++) {
	   /*JARKKO-damage*/
            iniat[iniind(nlayer,namorph,ninistate,n)].vx-=vxsum;
            iniat[iniind(nlayer,namorph,ninistate,n)].vy-=vysum;
            iniat[iniind(nlayer,namorph,ninistate,n)].vz-=vzsum;
         }
         MAXDEBUG(vxsum); MAXDEBUG(vysum); MAXDEBUG(vzsum);
      } /*initstate-looppi*/
     } /*amorph-looppi*/
    } /*layer-looppi*/
}

copyinilayer(iniat,at,atex,atflag,nlayer,ninistate,namorph,nat)
struct iniat *iniat;
struct atom at[];
struct atomex atex[];
struct atomflag atflag[];
int nlayer,ninistate,namorph,nat;
{
   int i;
   
   for (i=0;i<nat;i++) {
     /*JARKKO-damage*/
      at[i].x=iniat[iniind(nlayer,namorph,ninistate,i)].x;
      at[i].y=iniat[iniind(nlayer,namorph,ninistate,i)].y;
      at[i].z=iniat[iniind(nlayer,namorph,ninistate,i)].z;
      at[i].vx=iniat[iniind(nlayer,namorph,ninistate,i)].vx;
      at[i].vy=iniat[iniind(nlayer,namorph,ninistate,i)].vy;
      at[i].vz=iniat[iniind(nlayer,namorph,ninistate,i)].vz;
      at[i].ax=at[i].ay=at[i].az=0.0;
      at[i].a1x=at[i].a1y=at[i].a1z=0.0;
      at[i].a2x=at[i].a2y=at[i].a2z=0.0;
      at[i].type=iniat[iniind(nlayer,namorph,ninistate,i)].type;
      
      atex[i].pot=0.0;
      atflag[i].moving=True;
      atflag[i].interact=True;
      atflag[i].recoil=False;
      atflag[i].energetic=False;

   }
}

setinilayer(iniat,at,nlayer,ninistate,namorph,nat)
struct iniat *iniat;
struct atom at[];
int nlayer,ninistate,namorph,nat;
{
   int i;
   
   for (i=0;i<nat;i++) {
     /*JARKKO-damage*/
      iniat[iniind(nlayer,namorph,ninistate,i)].x=at[i].x;
      iniat[iniind(nlayer,namorph,ninistate,i)].y=at[i].y;
      iniat[iniind(nlayer,namorph,ninistate,i)].z=at[i].z;
      iniat[iniind(nlayer,namorph,ninistate,i)].vx=at[i].vx;
      iniat[iniind(nlayer,namorph,ninistate,i)].vy=at[i].vy;
      iniat[iniind(nlayer,namorph,ninistate,i)].vz=at[i].vz;
      iniat[iniind(nlayer,namorph,ninistate,i)].type=at[i].type;
      
   }
}




















