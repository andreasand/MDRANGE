#include "common.h"

initrec(file,reccalc,recresult,poly,recoil,unit,elstop,type)  

/* 
   Initialize entire recoil event calculation 
*/
struct file *file;
struct reccalc *reccalc;
struct recresult *recresult;
struct poly *poly;
struct recoil *recoil;
struct unit *unit;
struct elstop *elstop;
struct type *type;
{
   int i,j;
   real x;
   real getelstop();

   static char *str[] = {
      "unknown", 
      "default, surface at 0.0",
      "infinite negative backing"
   };

   for (i=0;i<7;i++) recresult->stopway[i]=0;
   poly->ndomain=poly->ndomainstot=0;
   poly->prevvel=0.0;
   poly->recpos.x=poly->recpos.y=poly->recpos.z=0.0;
   poly->range=0.0;
   if (reccalc->poly) {
      fprintf(file->fpo,"Polycrystalline mode selected\n");
      reccalc->type=2;
   }

   if (reccalc->type<1 || reccalc->type>2) 
      fprintf(file->fpdb,"Warning !!: illegal reccalc->type value %d\n",reccalc->type);

   if (reccalc->type != 1) 
      fprintf(file->fpo,"Recoil calc. type %d: %s\n\n",reccalc->type,str[reccalc->type]);

/* Print out recoil energy and velocity */

   recoil->attype=reccalc->Atype;
   reccalc->v0=sqrtf(2.0*reccalc->E0);
   fprintf(file->fpo,"Recoil energy %g eV and velocity %g m/s\n",
	   reccalc->E0,reccalc->v0*unit->vel[recoil->attype]);

   /*  check for implantations in Si JARKKO*/
   if(elstop->penr)
   {
   x=100*0.001*(reccalc->E0/type[recoil->attype].m)/(0.974*25);
   fprintf(file->fpo,"Recoil velocity is %g percent of V_Fermi of Si\n",x);
   if(x>50) fprintf(file->fpo,"WARNING! velocity may be too much for the PENR elstop theory!\n");
   }
   

/* Calculate the depen table step using a very simple method */

   recresult->dependz=reccalc->E0/1000/sqrt(type[recoil->attype].m);
   x=pow((double)10,floor(log10(recresult->dependz)));        
   recresult->dependz=x*floor((recresult->dependz+x/2.0)/x);      
   DEBUG(recresult->dependz);

   for (i=0;i<DEPENSIZE;i++) {
      recresult->depenel[i]=0.0;
      recresult->depennucl[i]=0.0;
   }
   recresult->sumFD=recresult->sumFDe=recresult->sumFDn=recresult->sumFDf=0.0; /*JARKKO*/

   recresult->drecspec=5.0;
   for (j=0;j<MAXTYPE;j++) {
      recresult->rectot[j]=recresult->rec0[j]=0.0;
      recresult->irectot[j]=recresult->irec0[j]=0;
      for (i=0;i<RECSPECSIZE;i++) {
	 recresult->recspec[j][i]=0.0;
	 recresult->irecspec[j][i]=0;
	 recresult->irecspec_Theta[j][i]=0; /*(JARKKO)*/
      }
   }
   /*JARKKO-damage*/
   for(i=0; i<ZBOXES; i++)
     recresult->depennucl_box[i]=0.0;

   for(i=0; i<DEPENSIZE; i++)
     {
     elstop->dens[i]=elstop->Sespc[i]=0.0;
     }

/* Initialize stopping table variables */

   recresult->stopvmax=reccalc->v0*1.1;  /* Leave some size above v0 */
   recresult->stopdv=recresult->stopvmax/STOPSIZE;
   if (file->printstop) {
      fprintf(file->fpo,"Stopping tables max velocity and velocity step %g %g\n",
	      recresult->stopvmax*unit->vel[recoil->attype],
	      recresult->stopdv*unit->vel[recoil->attype]);
      for (i=0;i<STOPSIZE;i++) {
	 recresult->elstop[i]=0.0;
	 recresult->nuclstop[i]=0.0;
	 recresult->ielstop[i]=0;
	 recresult->inuclstop[i]=0;
      }
   }

   /*JARKKO, this might be wrong, due to the layers-system for elstop ??? used layer 0 in getelstop*/
   if (file->printavvel) {
      /* 
	 Approximate the maximum t by using the limiting case
	 that the stopping is purely electronic with the S value
	 of the maximum velocity. 
      */
      recresult->avveltmax=2.5*reccalc->E0/(reccalc->v0*getelstop(elstop,unit,recoil,reccalc->v0/6.0,0));
      recresult->avveldt=recresult->avveltmax/AVVELSIZE;
      fprintf(file->fpo,"Average velocity tables maximum time and time step (fs) %g %g\n",
	      recresult->avveltmax*unit->time[recoil->attype]/1e-15,
	      recresult->avveldt*unit->time[recoil->attype]/1e-15);
      for (i=0; i<AVVELSIZE;i++) {
	 recresult->avvel[i]=0.0;
	 recresult->iavvel[i]=0;
	 recresult->curavvel[i]=0.0;
	 recresult->icuravvel[i]=0;
      }
   }

   recresult->nrange=0;

}

startrec(iniat,at,atex,atflag,layers,nlayer,ninistate,box,potcrit,time,
	 reccalc,recresult,poly,elstop,namorph)
/*
   Set some variables for start of a new recoil event
*/
struct iniat *iniat;
struct atom *at;
struct atomex *atex;
struct atomflag *atflag;
struct layers *layers;
int *nlayer,*ninistate,*namorph;
struct box *box;
struct potcrit *potcrit;
struct time *time;
struct reccalc *reccalc;
struct recresult *recresult;
struct poly *poly;
struct elstop *elstop;
{
   real unirand();
   int i;

   if (reccalc->poly) {
      getpolysize(poly,reccalc);
      poly->boundreached=False;
      if (poly->ndomain==0) {
	 poly->recpos.x=poly->recpos.y=poly->recpos.z=0.0;
	 poly->range=0.0;
      }
   }
   *nlayer=0;
   for (i=0;i<layers->Nlayers;i++) {
      layers->layer[i].crossed=False;
      layers->layer[i].selected=False;
   }
   *ninistate=unirand(0)*layers->Ninistate;
   if (*ninistate>=layers->Ninistate) *ninistate=layers->Ninistate-1;
   MAXDEBUGSRR("Getting inistate coords from (layer,inistate)",*nlayer,*ninistate);
   MAXDEBUGSRR("Getting inistate coords from (inistate,am.state)",*ninistate,*namorph);  /*JARKKO-damage*/
   copyinilayer(iniat,at,atex,atflag,*nlayer,*ninistate,*namorph,layers->layer[*nlayer].natoms[*namorph]); 
  
      /*JARKKO-damage */

   box->periodic=False;
   box->movement.x=0.0;
   box->movement.y=0.0;
   box->movement.z=0.0;

   potcrit->rm=potcrit->Rmrec;
   potcrit->r0=potcrit->R0rec;

   elstop->Se=0.0;elstop->totalfirsov=0.0;  /*JARKKO*/
   recresult->surfFDn=0.0;

   for (i=0;i<STOPSIZE;i++) {
      recresult->curelstop[i]=0.0;
      recresult->curnuclstop[i]=0.0;
      recresult->pcurnuclstop[i]=0.0; /*JARKKO*/
      recresult->icurelstop[i]=0;
      recresult->icurnuclstop[i]=0;
      recresult->picurnuclstop[i]=0; /*JARKKO*/
      recresult->TstopN[i]=0.0;
      recresult->TstopE[i]=0.0;
      recresult->iTstop[i]=0;
   }

   for (i=0; i<AVVELSIZE;i++) {
      recresult->curavvel[i]=0.0;
      recresult->icuravvel[i]=0;
   }

    recresult->stopN=recresult->stopE=recresult->stopv=recresult->stopR=0.0;
   recresult->istop=0;
  
}

createrecatom(reccalc,recoil,maxima,nat,at,atflag,poly,num)
/*
   Create the recoil atom 
*/
struct reccalc *reccalc;
struct recoil *recoil;
struct maxima *maxima;
int *nat;
struct atom *at;
struct atomflag *atflag;
struct poly *poly;
int num;
{
   int i,seed1,j,k,l;
   real fii,theta,t;
   static real pfii,ptheta; /* Previous polycr. domain alignment angle */
   real unirand();
   real newx,newy,newz,newvx,newvy,newvz;
   int gridsize;

   (*nat)++;
   i=*nat-1;

   at[i].ax=at[i].ay=at[i].az=0.0;
   at[i].a1x=at[i].a1y=at[i].a1z=0.0;
   at[i].a2x=at[i].a2y=at[i].a2z=0.0;
   at[i].type=reccalc->Atype;

   atflag[i].moving=atflag[i].recoil=atflag[i].energetic=True;

   recoil->number=i; recoil->attype=at[i].type;
   recoil->splitnum=0; /*reed*/

   if (debug) unirand(-1);
   seed1=0;
   if (reccalc->Forceseed!=0) seed1=reccalc->Forceseed;

   recoil->ddx=recoil->ddy=recoil->ddz=0.0;


   /*old random*/
   
   recoil->x=recoil->x0=at[i].x=reccalc->Startmin.x+unirand(seed1)*
                                   (reccalc->Startmax.x-reccalc->Startmin.x);
   recoil->y=recoil->y0=at[i].y=reccalc->Startmin.y+unirand(0)*
                                   (reccalc->Startmax.y-reccalc->Startmin.y);
   recoil->z=recoil->z0=at[i].z=reccalc->Startmin.z+unirand(0)*
                                   (reccalc->Startmax.z-reccalc->Startmin.z);
   
  

/*
   Select recoil fii and theta. Take care with possible polycrystallinity modes !
*/
   if (reccalc->poly && poly->ndomain != 0 && reccalc->Polydeltaangle > 0.0) {

      /* Polycrystalline mode theta,phi  selection */

      recoil->fii = fii = pfii+(2*unirand(0)-1.0)*reccalc->Polydeltaangle;
      theta = ptheta+(2*unirand(0)-1.0)*reccalc->Polydeltaangle;
      /* Now beware: |theta| must not be greater than 45 degrees. But ptheta can be
	 anything ! So, we have to transform theta by adding 90 degrees to it ! */
      while (oabs(theta) > 3.1415927/4.0+0.001) { /* Some tolerance to avoid round-off trouble ! */
	 theta-=sign(theta)*3.1415927/2.0;
      }
      recoil->theta=theta;
      DEBUGSRRRR("Polycr. angles",pfii,fii,ptheta,theta);
   }
   else {

      /* Ordinary theta, phi selection */

      if (reccalc->Fiimax==0.0 && reccalc->Thetamax==0.0) {
	 recoil->fii  = fii  = reccalc->Fii0;
	 recoil->theta= theta= reccalc->Theta0;
      }
      else {
	 recoil->fii=fii=reccalc->Fii0+unirand(0)*(reccalc->Fiimax-reccalc->Fii0);
	 if (reccalc->Theta0 == 0.0) {
	    t=reccalc->Thetamax;
	    recoil->theta=theta=acos(cos(t)+unirand(0)*(1.0-cos(t)));
	 }
	 else {
	    recoil->theta=theta=reccalc->Theta0+unirand(0)*(reccalc->Thetamax-reccalc->Theta0);
	 }
	 if (unirand(0) < 0.5) recoil->theta=theta=-theta;
      }
   }
   pfii=fii;       
   ptheta=theta;

/*
   Initial pos. and direction selected, now set velocities
*/

   MAXDEBUGSRRR("Recoil start",recoil->x0,recoil->y0,recoil->z0);
   MAXDEBUGSRR("Recoil start",recoil->fii,recoil->theta);

   recoil->vel=reccalc->v0=sqrtf(2.0*reccalc->E0);;
   if (reccalc->poly && poly->ndomain != 0) {
      recoil->vel=poly->prevvel;
   }

   maxima->vel=recoil->vel;
   maxima->vel_attype=recoil->attype;

   recoil->energy=0.5*recoil->vel*recoil->vel;
   recoil->eprev=recoil->energy; recoil->vprev=recoil->vel;
   DEBUGSRR("Giving new recoil",recoil->vel,recoil->energy);
   recoil->vel0=recoil->vel;
   at[i].vx=recoil->vx0=recoil->vel*sin(theta)*cos(fii);
   at[i].vy=recoil->vy0=recoil->vel*sin(theta)*sin(fii);
   at[i].vz=recoil->vz0=recoil->vel*cos(theta);

   
   MAXDEBUGSRRRR("Recoil start", recoil->vx0,recoil->vy0,recoil->vz0,recoil->vel0);

/*
   Set some parameters at start of calculation
*/
   recoil->acc=0.0;
   maxima->acc=recoil->acc;
   maxima->acc_attype=recoil->attype;

   recoil->touched=False; /*JARKKO surface debug*/
   recoil->beeninside=False;
   recoil->realbeeninside=False; /*JARKKO-damage*/
   recoil->turnedaround=False;
   recoil->incollision=False;
   recoil->range=recoil->r_path=recoil->r_chord=recoil->r_proj=0.0;
   
   recoil->n++;
   recoil->xp=recoil->x; recoil->yp=recoil->y; recoil->zp=recoil->z;
   recoil->sx=recoil->x; recoil->sy=recoil->y; recoil->sz=recoil->z;

   recoil->r_proj=(recoil->vx0*recoil->sx+recoil->vy0*recoil->sy+
		   recoil->vz0*recoil->sz)/recoil->vel0;
   
}

getrecoilprop(recoil,reccalc,at,box,poly)
struct recoil *recoil;
struct reccalc *reccalc;
struct atom *at;
struct box *box;
struct poly *poly;
{
   int i;
   real r;

   i=recoil->number;

   recoil->x=at[i].x;
   recoil->y=at[i].y;
   recoil->z=at[i].z;
   recoil->vx=at[i].vx;
   recoil->vy=at[i].vy;
   recoil->vz=at[i].vz;
   recoil->sx=at[i].x+box->movement.x;  recoil->px=recoil->sx+poly->recpos.x;
   recoil->sy=at[i].y+box->movement.y;  recoil->py=recoil->sy+poly->recpos.y;
   recoil->sz=at[i].z+box->movement.z;  recoil->pz=recoil->sz+poly->recpos.z;

   recoil->dx=recoil->sx-recoil->xp; 
   recoil->dy=recoil->sy-recoil->yp; 
   recoil->dz=recoil->sz-recoil->zp; 
   recoil->dr=sqrtf(recoil->dx*recoil->dx+recoil->dy*recoil->dy+recoil->dz*recoil->dz);

   recoil->xp=recoil->sx;
   recoil->yp=recoil->sy;
   recoil->zp=recoil->sz;

   recoil->eprev=recoil->energy;
   recoil->vprev=recoil->vel;
   recoil->vel=sqrtf(at[i].vx*at[i].vx+at[i].vy*at[i].vy+at[i].vz*at[i].vz);
   recoil->acc=sqrtf(at[i].ax*at[i].ax+at[i].ay*at[i].ay+at[i].az*at[i].az);
   recoil->energy=0.5*recoil->vel*recoil->vel;
   
   recoil->r_path+=recoil->dr;
   recoil->r_chord=sqrtf((recoil->x0-recoil->sx)*(recoil->x0-recoil->sx)+
			 (recoil->y0-recoil->sy)*(recoil->y0-recoil->sy)+
			 (recoil->z0-recoil->sz)*(recoil->z0-recoil->sz));
   r=recoil->r_proj;
   recoil->r_proj=(recoil->vx0*recoil->sx+recoil->vy0*recoil->sy+
		   recoil->vz0*recoil->sz)/recoil->vel0;
   recoil->drp=recoil->r_proj-r;

   if (reccalc->Trange == 0) r=recoil->sz;
   if (reccalc->Trange == 1) r=recoil->r_proj;
   if (reccalc->Trange == 2) r=recoil->r_chord;
   if (reccalc->Trange == 3) r=recoil->r_path;
   recoil->range=r;

   if (at[i].vz < 0) recoil->turnedaround=True;
   MAXDEBUGSRRRR("getrecprop:",r,recoil->vel,recoil->acc,recoil->dr);

}



getmaximumprop(maxima,at,nat)
/*
    Get maximum properties for time step calculation
    We need vel, acc and attype
*/
struct maxima *maxima;
struct atom *at;
int nat;
{
   int i;
   real vel,acc;

   maxima->vel=0; maxima->acc=0;
   for (i=0;i<nat;i++) {
     
     vel=(at[i].vx*at[i].vx+at[i].vy*at[i].vy+at[i].vz*at[i].vz);
     if (vel>maxima->vel) { maxima->vel=vel; maxima->vel_attype=at[i].type; }
     
     acc=(at[i].ax*at[i].ax+at[i].ay*at[i].ay+at[i].az*at[i].az);
     if (acc>maxima->acc) { maxima->acc=acc; maxima->acc_attype=at[i].type; }
   }
   maxima->vel=sqrtf(maxima->vel); /*(JARKKO)*/
   maxima->acc=sqrtf(maxima->acc); /*(JARKKO)*/

   MAXDEBUGSRR("getmaxprop:",maxima->vel,maxima->vel_attype);
   MAXDEBUGSRR("getmaxprop:",maxima->acc,maxima->acc_attype);

}


calcmisc(recoil,elstop,reccalc,recresult,poly,file,firstrecstep,t,physical,nreed)
/*
  Integrating the deposited energy distribution
  in the end of the calculation
  yields a somewhat (about 4 %) too small value. This probably is because
  the dividing with dependz is a bit wrong at the surface and the
  end of the calculation, and could be corrected by making dependz smaller.

  In poly calculation, poly->range must be taken in account !
*/
struct recoil *recoil;
struct elstop *elstop;
struct reccalc *reccalc;
struct recresult *recresult;
struct poly *poly;
struct file *file;
struct physical *physical;
int nreed;
logical firstrecstep;
real t;
{
   double dr,r,FD,FDn,FDe,FDf;  
   double S,Sn,Se;
   int i;

   static double Einstartofcollision,Sesum,drsum;
   static int nstepsincollision;
   static logical incollisionlasttime;
   static logical warnavveloutside;

   if (firstrecstep) {
      Sesum=0.0;
      drsum=0.0;
      nstepsincollision=0;
      incollisionlasttime=False;
      Einstartofcollision=-1.0;   /* Error flag */
      warnavveloutside=False;
   }

   dr=recoil->dr; /*dz?? */
   if (reccalc->Trange == 1) dr=recoil->drp;
   /* if (dr == 0.0) return; */

   FD=(recoil->eprev-recoil->energy);
   
   FDe=elstop->Se*recoil->dr; /* * pow(0.5,recoil->splitnum); virtual parts?*/
   FDf=elstop->totalfirsov*recoil->dr;/* * pow(0.5,recoil->splitnum); virtual parts?*/
   FDn=FD-FDe-FDf; 

   
   
   if (file->printdepen) {
     if(reccalc->reed)
       {
	 FD*=pow(0.5,recoil->splitnum);
	 FDn*=pow(0.5,recoil->splitnum);
	 FDe*=pow(0.5,recoil->splitnum);
	 FDf*=pow(0.5,recoil->splitnum);
       }
      r=recoil->sz;
      if (reccalc->Trange == 1) r=recoil->r_proj;
      if (reccalc->poly) r+=poly->range;

      i=(int) (r/recresult->dependz+0.5);
      if (i>DEPENSIZE) i=DEPENSIZE;
      if (i<0) i=0;
        recresult->depennucl[i]+=FDn;
        recresult->depenel[i]+=FDe;
        recresult->sumFD+=FD;
        recresult->sumFDn+=FDn;
        recresult->sumFDe+=FDe;	
	recresult->sumFDf+=FDf; /*JARKKO firsov*/
   }

   MAXDEBUGSRRRRR("calcdepen.c",dr,i,FDn,elstop->Se,recoil->energy);
  
   if (file->printstop && nreed==0) {
      /*
	 Note that the calculation of the nuclear stopping power is a real
	 problem ! Since the time step is always calculated from the previous
	 force, there exists an asymmetry in the length of delta r in collisions
	 which makes the calculation of Sn difficult !? I believe this can be 
	 remedied by calculating the stopping over an entire collision

	 Also, one must remember not to include stopping of discarded particles !

	 But the _curve_ of S is not just dE/dz, but actually (dE/dz)dv ! This 
	 must be so since the intervals are over v !?!? Bullshit, I think.
      */
      i=(int) recoil->vel/recresult->stopdv;
      Se=Sn=S=0.0;
      if (i>=0 && i<STOPSIZE) {
	 if (dr != 0.0) {
	   S=(recoil->eprev-recoil->energy)/dr;/**(recoil->vprev-recoil->vel);*/
	    Se=elstop->Se;
	    Sn=(S-Se);
	    
	    /*printf("S:%lg Sn=%lg Se=%lg\n",S,Sn,Se);*/
	    /*JARKKO, only positive values of Sn are calculated also */
	    if(Sn>0){recresult->pcurnuclstop[i]+=Sn;
	             recresult->picurnuclstop[i]++;}
	    recresult->curnuclstop[i]+=Sn;
	    recresult->icurnuclstop[i]++;
	    recresult->curelstop[i]+=Se;
	    recresult->icurelstop[i]++;

	    /*kusee*/
	    
	    recresult->stopN+=Sn;
	    recresult->stopE+=Se;
	    recresult->stopv+=recoil->vel;
	    recresult->istop++;
	    recresult->stopR+=recoil->dr;
	    
	    /*printf("kulku: %lg(%d) sn:%lg se:%lg\n",recresult->stopR,recresult->istop,recresult->stopN,recresult->stopE);*/

	    if(recresult->stopR>reccalc->stopDR)
	      {
	      i=(int) (recresult->stopv/recresult->istop) /recresult->stopdv;  /*average velocity in the dr-slab.*/
	      recresult->TstopN[i]+=recresult->stopN/recresult->istop;
	      recresult->TstopE[i]+=recresult->stopE/recresult->istop;
	      recresult->iTstop[i]++;	      
	      recresult->stopN=recresult->stopE=recresult->stopv=recresult->stopR=0.0;
	      recresult->istop=0;	    	      
	      }

	 }

	 MAXDEBUGSRRR("stopping",i,Sn,Se);
	 MAXDEBUGSRRR("stopping",i,recresult->curnuclstop[i],recresult->curelstop[i]);
	 MAXDEBUGSRRR("stopping",i,recresult->icurnuclstop[i],recresult->icurelstop[i]);

      }

   }
   
   if (file->printavvel && nreed==0) {
      i=t/recresult->avveldt;
      if (i>=0 && i<AVVELSIZE) {
	/*if(reccalc->reed)
	  count=pow(0.5,recoil->splitnum);
	  else count=1;*/
	recresult->curavvel[i]+=recoil->vel;
	recresult->icuravvel[i]++;
      }
      else {
	 if (!warnavveloutside) {
	    fprintf(fpdebug,"calcmisc: too small borders for avvel calc.: %d %g\n",i,t);
	    warnavveloutside=True;
	 }
      }
      MAXDEBUGSRRRR("avvel",i,t,recoil->vel,recresult->avvel[i]);
   }

}


endrec(at,recoil,reccalc,recresult,poly,file,nat,nrec,recerror,physical,nreed) 
struct atom *at;
struct recoil *recoil;
struct reccalc *reccalc;
struct recresult *recresult;
struct poly *poly;
struct file *file;
int *nat,nrec;
logical recerror;
struct physical *physical;
int nreed;
{
   double Snsum,Sesum;
   real r,E;
   int i,ind,itype;
   /*JARKKO*/
   int numreal;
   real table[MAXRECOIL];
/*
   Stopway holds the way the recoil has stopped:
   stopway = 0:  Energy below threshold reccalc->Emin 
           = 1:  Escaped (sputtered)
	   = 2:  Error !! (something completely insane)
	   = 3:  Range outside array for printing in printresults;
	   = 4:  as 0, but negative, if negative is allowed (reccalc->type==2)
	   = 6:  Domain boundary reached in polycrystalline calc.
	   = 7:  Zmax reached in calculation

   Valid stopways for the range being stored in recresult are 0 and 4.

*/
   int stopway;

   (*nat)--;

   stopway=0;

   if (reccalc->Trange == 0) r=recoil->sz;
   if (reccalc->Trange == 1) r=recoil->r_proj;
   if (reccalc->Trange == 2) r=recoil->r_chord;
   if (reccalc->Trange == 3) r=recoil->r_path;
   recoil->range=r;

   if (recoil->range < 2.0*reccalc->Startmin.z) stopway=1;
   if (reccalc->type==2 && recoil->range < 0.0) stopway=4;
   if (recoil->range>reccalc->Zmax) stopway=7;

   if (r < -1e20 || r > 1e20) stopway=2;  /* Attempt to localize NaN:s */
   if (r != r) stopway=2;  /* IEEE standard defines that x != x if x is NaN */
   if (recerror) stopway=2;

   if (file->printrecend) {
     fprintf(file->fpo,"Recoil end pos %g %g %g E %g r %g\n",
	     recoil->sx,recoil->sy,recoil->sz,recoil->energy,recoil->range);
   }


   if (reccalc->poly) { /* End of polycrystalline domain */
      poly->recpos.x+=recoil->sx; 
      poly->recpos.y+=recoil->sy; 
      poly->recpos.z+=recoil->sz;
      poly->range+=recoil->r_proj;
      if (poly->boundreached && !recerror) {
	 stopway=6;
      }
      else {   /* End of polycrystalline range calc. */
	 poly->ndomain=0;
	 recoil->range=poly->range;
	 if (recoil->range<0 && !recerror) stopway=4;
	 else if (!recerror) stopway=0;
      }
      DEBUGSRRR("endrec.c poly end",poly->ndomain,poly->range,stopway);
      MAXDEBUGSRRR("endrec.c poly end",poly->recpos.x,poly->recpos.y,poly->recpos.z);
   }

   MAXDEBUGSRRRRR("endrec.c",stopway,r,recoil->sx,recoil->sy,recoil->sz);
   recresult->stopway[stopway]++;
   if (stopway == 0 || stopway == 4) {
      if (recresult->nrange >= MAXRECOIL) {
	printf("mdh ERROR: MAXRECOIL overflow. Exiting.\n");
	exit(0);
      }

      
      if(physical->sputtering) { /*JARKKO*/
	numreal=0;
	for(i=0; i<recresult->nrange; i++)
	  {
	  recresult->ranges[i]-=physical->deltasput_surf;  /*edellisista vahennetaan juuri irronnut aines*/
	  if(recresult->ranges[i]>=0) {table[numreal]=recresult->ranges[i];numreal++;}
	  }
      for(i=0; i<numreal; i++)
	recresult->ranges[i]=table[i];
      recresult->nrange=numreal;   /*pinnan ulkopuolelta poistetaan..*/
      }
      recresult->ranges[recresult->nrange]=recoil->range- physical->sput_surf; /*JARKKO erosion,range is counted from surface*/
      recresult->weights[recresult->nrange]=pow(0.5,recoil->splitnum); /*sputtering indexing goes wrong.*/
      recresult->nrange++;

      if (file->printstop) {
	 Snsum=Sesum=0.0;
	 for (i=0;i<STOPSIZE;i++) {
	    recresult->nuclstop[i]+=recresult->curnuclstop[i];
	    recresult->inuclstop[i]+=recresult->icurnuclstop[i];
	    recresult->pnuclstop[i]+=recresult->pcurnuclstop[i];
	    recresult->pinuclstop[i]+=recresult->picurnuclstop[i];
	    recresult->elstop[i]+=recresult->curelstop[i];
	    recresult->ielstop[i]+=recresult->icurelstop[i];
	    if (recresult->icurnuclstop[i] != 0) {
	       Snsum+=recresult->curnuclstop[i]/recresult->icurnuclstop[i];
	       Sesum+=recresult->curelstop[i]/recresult->icurelstop[i];
	    }
	 }
	 DEBUGSRR("Stoppings for this event",Snsum,Sesum);
      }

      if (file->printavvel) {
	 for(i=0;i<AVVELSIZE;i++) {
	    if (recresult->icuravvel[i] != 0) {
	       recresult->avvel[i]+=recresult->curavvel[i]/recresult->icuravvel[i];
	       recresult->iavvel[i]+=recresult->icuravvel[i];
	    }
	 }
      }
   }

   if (stopway==1) recresult->sumFDsput+=recoil->energy; /*note that, in reality ions have to overcome the binding energy barrier to get sputtered!*/
   if (stopway==7) recresult->sumFDtrans+=recoil->energy;

   if (file->recoilspec && nreed==0) {
      /* Add up energy of now final atoms into recoil spectrum */
      for (i=0;i< *nat;i++) {
	 if (i==recoil->number) continue;
	 E=0.5*(at[i].vx*at[i].vx+at[i].vy*at[i].vy+at[i].vz*at[i].vz);
	 itype=at[i].type;
	 MAXDEBUGSRR("recspec",i,E);
	 recresult->rectot[itype]+=E;
	 recresult->irectot[itype]++;
	 
	 if (E>exp(-0.5/recresult->drecspec)) {
	    ind=(int)(log(E)*recresult->drecspec+0.5);
	    if (ind > RECSPECSIZE) {
	      printf("mdh endrec ERROR: RECSPECSIZE overflow %d\n",ind);
	      exit(0);
	    }
	    recresult->recspec[itype][ind]+=E;
	    recresult->irecspec[itype][ind]++;
	 }
	 if (E<1.0) {
	    recresult->rec0[itype]+=E;
	    recresult->irec0[itype]++;
	 }
      }
   }


   
}
printresults(file,reccalc,recresult,unit,physical,box,elstop,nrec)
struct file *file;
struct reccalc *reccalc;
struct recresult *recresult;
struct unit *unit;
struct physical *physical;
struct box *box;
struct elstop *elstop;
int nrec;
{
   int i,j,n,jmax,itype,ityperecmax;
   real min,minr,max,maxr,sum,mean,deltar,sign;
   real dmean,sigmar,straggle;
   real r,d,x,Erectot,Erec0,dE;
   real rarray[RARRAYSIZE];
   double Snsum,Sesum;
   double avvelsum;
   char buf[120];
   real a,b,area,count; /*JARKKO*/

   if (nrec<2) return;

   fprintf(file->fpo,"Output of %d recoil events at ",nrec); eltime();


   for (i=0;i<RARRAYSIZE;i++) rarray[i]=0;

   min=1e20; max=-1e20; sum=0.0;
   mean=0.0; sigmar=0.0; straggle=0.0; count=0.0;
   for(i=0;i<recresult->nrange;i++) {
      r=recresult->ranges[i];
      if (min > r) min=r;
      if (max < r) max=r;
      if(reccalc->reed) {sum+=r*recresult->weights[i];count+=recresult->weights[i];}
      else {sum+=r;count++;}
   }
   if (recresult->nrange>0) mean=sum/count;
   DEBUGSRRRR("range vars",recresult->nrange,min,max,mean);

   sign=1;
   deltar=MIN(mean-min,max-mean)/15;      /* Get first value of printing interval */
   if (deltar <= 1e-10 || deltar >= 1e10) { 
      fprintf(fpdebug,"Warning: deltar %g min %g mean %g max %g\n",deltar,min,mean,max);
      fprintf(fpdebug,"Setting deltar to 1.0\n");
      deltar=1.0;
   }
   x=pow((double)10,floor(log10(deltar)));        /* Get power of 10 less than deltar */
   deltar=x*floor((deltar+x/2.0)/x);      /* Round interval to one digit precision */
   minr=deltar*floor(min/deltar);         /* Get minimum */

   jmax=-1;
   recresult->stopway[3]=0;
   dmean=0.0;
   jmax=1;
   for (i=0;i<recresult->nrange;i++) {      
     if(reccalc->reed) {dmean+=(recresult->ranges[i]*recresult->weights[i]-mean)*(recresult->ranges[i]*recresult->weights[i]-mean);}
      else dmean+=(recresult->ranges[i]-mean)*(recresult->ranges[i]-mean);

      r=recresult->ranges[i]+deltar/2.0-minr;
      j=(int) floor(r/deltar);
      if (j>=RARRAYSIZE || j<0) { 
	 recresult->stopway[3]++;
	 if (debug) { 
	    fprintf(fpdebug,"printresults.c Warning: range outside rarray r %g j %d\n",
		    recresult->ranges[i],j);
	    fprintf(fpdebug,"min %g mean %g max %g deltar %g minr %g\n",
		    min,mean,max,deltar,minr);
	 }
	 continue;  /* Don't modify jmax or rarray if out of range ! */
      }
      if (j>jmax) jmax=j;
      rarray[j]+=recresult->weights[i];     
   }

   if (count>1) {  /*oli 0, mutta kaatumisen estamiseksi 1 (JARKKO)*/
     straggle=sqrtf(dmean/count);
     sigmar=sqrtf(dmean/count/(count-1));
   }

   fprintf(file->fpo,"Ways of stopping:");
   for (i=0;i<=7;i++) fprintf(file->fpo," %d",recresult->stopway[i]);
   fprintf(file->fpo,"\n");

   fprintf(file->fpo,"Mean range %g +- %g straggle %g\n",mean,sigmar,straggle);

   /* Range distribution output */

   file->fpra=fopen(file->rangef,"w");
   for (j=0;j<=jmax;j++) fprintf(file->fpra,"%g %g\n",minr+j*deltar,rarray[j]);
   fclose(file->fpra);

   /*real distribution according to the dose (JARKKO) */
   
   /*calculate the area*/
   area=0;
   for(j=0; j<=jmax; j++)
      {
      if(j==0) a=0;
      else a=minr+(j-1)*deltar+deltar/2;
      b=minr+(j+1)*deltar-deltar/2;
      area+=(b-a)*rarray[j];
      }
   if(area==0) {area=1;printf("No range data at all!! realrange.out is impossible to make.\n");}
   else printf("Area for range.out is %g\n",area);
   file->fpra=fopen("realrange.out","w");
   for (j=0;j<=jmax;j++) fprintf(file->fpra,"%g %g\n",minr+j*deltar,1e8*physical->Dose*(float)rarray[j]/area);  /*1e8 comes from unit transv. realrange -> [1/cm^3] */
   fclose(file->fpra);

   /* Deposited energies output */

   if (file->printdepen) {
      fprintf(file->fpo,"Total deposited energies/ion: FDn %g FDe %g FDf %g\n",
	      recresult->sumFDn/nrec,recresult->sumFDe/nrec,recresult->sumFDf/nrec);
      for(i=DEPENSIZE-1;i>1;i--) if (recresult->depennucl[i]!=0.0) break;
      DEBUGSR("Writing depen file",i);
      file->fpde=fopen(file->depenf,"w");
      for (j=0;j<=i;j++) {
         d=nrec*recresult->dependz;
         if (j==0) d/=2.0;   /* Surface at 0 ! */
	 fprintf(file->fpde,"%g %g %g\n",j*recresult->dependz,
		 recresult->depennucl[j]/d,recresult->depenel[j]/d);
      }
      fclose(file->fpde);
   }

   if(physical->Damage)
     {
     /*JARKKO-damage*/
       
      for(i=ZBOXES-1;i>0;i--) if(recresult->depennucl_box[i]!=0.0) break;
     
      file->fpde=fopen("boxdepen.out","w");
      Sesum=0.0; /*kokonaissumma*/
      for (j=0;j<=i;j++) {
	
	Sesum+=recresult->depennucl_box[j];
	fprintf(file->fpde,"%g %g \n",j*box->movelim.z,recresult->depennucl_box[j]*physical->Dosescaling);
      }/*Ong vs. eV/Ong^3*/

      /*fprintf(file->fpde,"%g\n",Sesum);*/
		
      fclose(file->fpde);
     
     }
   
   fprintf(file->fpo,"Energy/ion lost outside: sputtered %g transmitted %g eV/ion\n",
	   recresult->sumFDsput/nrec,recresult->sumFDtrans/nrec);

   /* Recoil spectrum output */

   if (file->recoilspec) {
      for (itype=0;itype<MAXTYPE;itype++) if (recresult->irectot[itype]!=0) ityperecmax=itype;
      for (itype=0;itype<MAXTYPE;itype++) {
	 if (recresult->irectot[itype]==0) continue;
	 /* Remember to subtract out the thermal energy from the total energy results */
	 Erectot=recresult->rectot[itype]-1.5*recresult->irectot[itype]*physical->Trec/unit->tempconv;
	 Erec0=recresult->rec0[itype]-1.5*recresult->irec0[itype]*physical->Trec/unit->tempconv;
      
	 fprintf(file->fpo,"Energy/ion lost to %.2f recoils of type %d: total %g < 1 eV %g\n",
		 1.0*(recresult->irectot[itype]-recresult->irec0[itype])/nrec,itype,
		 Erectot/nrec,Erec0/nrec);
      
	 for (i=RECSPECSIZE-2;i>1;i--) if (recresult->irecspec[itype][i]!=0) break;
	 DEBUGSRR("Writing recoil spectrum file",i,itype);
	 if (ityperecmax!=1) sprintf(buf,"%s.%d.out",file->recspecf,itype);
	 else sprintf(buf,"%s.out",file->recspecf,itype);
	 file->fprs=fopen(buf,"w");
	 for (j=0;j<=i+1;j++) {
	    d=exp(j/recresult->drecspec);
	    dE=exp((j+0.5)/recresult->drecspec)-exp((j-0.5)/recresult->drecspec);
	    fprintf(file->fprs,"%10.2f %12.6g %10.3f %10.4f %10.6g %9d\n",d,
		 1.0*recresult->irecspec[itype][j]/nrec/dE,dE,
		 1.0*recresult->irecspec[itype][j]/nrec,
		    recresult->recspec[itype][j]/nrec,
		    recresult->irecspec[itype][j]);
	 }
	 fclose(file->fprs);
	 /*(JARKKO)*/
	 for (i=RECSPECSIZE-2;i>1;i--) if (recresult->irecspec_Theta[itype][i]!=0) break;
	 DEBUGSRR("Writing angular recoil spectrum file",i,itype);
	 if (ityperecmax!=1) sprintf(buf,"%s.%d.out",file->recspec_Thetaf,itype);
	 else sprintf(buf,"%s.out",file->recspec_Thetaf,itype);
	 file->fprst=fopen(buf,"w");
	 fprintf(file->fprst,"Lowest energy:%i, Highest energy:%i \n",reccalc->ETmin,reccalc->ETmax);
	 for(j=0;j<=i+1;j++) 
	    {
	    fprintf(file->fprst,"%i %i\n",j,recresult->irecspec_Theta[itype][j]);
	    }
	 fclose(file->fprst);
	 /*(JARKKO)*/
      }
   }


   /* Stopping output */

   if (file->printstop) {
      DEBUGS("Writing stopping tables file\n");
      file->fpst=fopen(file->stopf,"w");

      fprintf(file->fpst,"%g %g %g\n",0.0,0.0,0.0);
      Snsum=0.0; Sesum=0.0;
      for(i=0;i<STOPSIZE;i++) {
	 if (recresult->inuclstop[i]>0) {
	    Snsum+=recresult->nuclstop[i]/recresult->inuclstop[i];
	    Sesum+=recresult->elstop[i]/recresult->ielstop[i];
	    if(recresult->iTstop[i]==0) recresult->iTstop[i]=1;
	    fprintf(file->fpst,"%g %g %g %g %g %g\n",
		    (i+0.5)*recresult->stopdv*unit->vel[reccalc->Atype],
		    recresult->nuclstop[i]/recresult->inuclstop[i],
		    recresult->elstop[i]/recresult->ielstop[i],
		    recresult->pnuclstop[i]/recresult->pinuclstop[i],
		    recresult->TstopN[i]/recresult->iTstop[i],
		    recresult->TstopE[i]/recresult->iTstop[i]);
	 }
	 else {
	    fprintf(file->fpst,"%g %g %g\n",(i+0.5)*
		    recresult->stopdv*unit->vel[reccalc->Atype],0.0,0.0);
	 }
      }
      fclose(file->fpst);
      fprintf(file->fpo,"Total stoppings: Sn %g Se %g\n",Snsum,Sesum);

   }

    if (file->printdens) {
      DEBUGS("Writing dens spectrum file\n");
      file->fpst=fopen("eldens_spec.out","w");
      for(i=0;i<DEPENSIZE;i++) 
	if (elstop->dens[i]>0 || elstop->Sespc[i]>0) 
	   fprintf(file->fpst,"%g %g %g %g\n",(float)i*(10.0/DEPENSIZE),(float)i*(10.0/DEPENSIZE)/0.1481848229
		   ,elstop->dens[i]/nrec,elstop->Sespc[i]/nrec);
      fclose(file->fpst);
    } /*density: e/a0^3, e/Ong^3, average count per ion, average Se per ion*/

   /* Average velocity output */

   if (file->printavvel) {
      DEBUGS("Writing average velocity file\n");
      avvelsum=0.0;
      file->fpav=fopen(file->avvelf,"w");
      for(i=0; i<AVVELSIZE;i++) {
	 if (recresult->iavvel[i] !=0 && recresult->nrange>0 ) {
	    fprintf(file->fpav,"%g %g %d\n",
		    (i+0.5)*recresult->avveldt*unit->time[reccalc->Atype]/1e-15,
		    recresult->avvel[i]/recresult->nrange*unit->vel[reccalc->Atype],
		    recresult->iavvel[i]*recresult->nrange);
	    avvelsum+=recresult->avvel[i]/recresult->nrange*unit->vel[reccalc->Atype]*
	       recresult->avveldt*unit->time[reccalc->Atype]/1e-15;
	 }
      }
      fclose(file->fpav);

      fprintf(file->fpo,"Integrated average velocity %g\n",avvelsum);
   }
   fprintf(file->fpo,"\n");

}
/*----------------------------------------------------------------------------*/
getestat(recoil,reccalc,file,printflag)
struct recoil *recoil;
struct reccalc *reccalc;
struct file *file;
logical printflag;
{
   static int *estat,*fiistat,*thetastat;
   static logical firsttime=True;
   int i,j;
   int eind(),find(),tind();
   real zstep,estep,fstep,tstep;
   real zmax,thetamax=180,fiimax=360;
   real fii,theta;
   real rad2deg;

   if (firsttime==True) {
      fprintf(fpdebug,"Warning: getestat doesn't work yet !!!\n");
      firsttime=False;
      estat=malloc(ESTATZSIZE*ESTATESIZE*sizeof(int));
      fiistat=malloc(ESTATZSIZE*ESTATFIISIZE*sizeof(int));
      thetastat=malloc(ESTATZSIZE*ESTATTHETASIZE*sizeof(int));
      for (i=0; i<ESTATZSIZE; i++) {
	 for (j=0;j<ESTATESIZE;j++) estat[eind(i,j)]=0;
	 for (j=0;j<ESTATFIISIZE;j++) fiistat[find(i,j)]=0;
	 for (j=0;j<ESTATTHETASIZE;j++) thetastat[tind(i,j)]=0;
      }
   }

   rad2deg=180.0/3.1415927;
   zmax=reccalc->Estatzmax;
   if (!printflag) {
      j=(int) (recoil->x/zmax*ESTATZSIZE);
      j=WITHIN(j,0,ESTATZSIZE-1);
      /* Calculate energy index */
      i=(int) (recoil->energy/reccalc->E0*ESTATESIZE);
      i=WITHIN(i,0,ESTATESIZE-1);
      estat[eind(i,j)]++;
      /* Calculate fii index (v dot x = v x cos(fii) */
      fii=acos(recoil->vx/recoil->vel);
      i=(int) fii*rad2deg/180.0*ESTATFIISIZE;
      i=WITHIN(i,0,ESTATFIISIZE-1);
      estat[eind(i,j)]++;
      /* Calculate theta index (v dot z = v z cos(theta) */
      theta=acos(recoil->vz/recoil->vel);
      i=(int) theta*rad2deg/180.0*ESTATTHETASIZE;
      i=WITHIN(i,0,ESTATTHETASIZE-1);
      estat[eind(i,j)]++;      
      
   }
   else {   /* Print out arrays */
      file->fpee=fopen(file->estatef,"w");
      file->fpef=fopen(file->estatff,"w");
      file->fpet=fopen(file->estattf,"w");
      zstep=reccalc->Estatzmax/ESTATZSIZE;
      for (i=0; i<ESTATZSIZE; i++) {
	 fprintf(file->fpee,"%g",i*zstep);
	 fprintf(file->fpef,"%g",i*zstep);
	 fprintf(file->fpet,"%g",i*zstep);
	 for (j=0;j<ESTATESIZE;j++) fprintf(file->fpee," %d",estat[eind(i,j)]);
	 for (j=0;j<ESTATFIISIZE;j++) fprintf(file->fpef," %d",fiistat[find(i,j)]);
	 for (j=0;j<ESTATTHETASIZE;j++) fprintf(file->fpet," %d",thetastat[tind(i,j)]);
	 fprintf(file->fpee,"\n");
	 fprintf(file->fpef,"\n"),
	 fprintf(file->fpet,"\n");
      }

   }
}
int eind(i,j)
int i,j;
{
   return(i*ESTATESIZE+j);
}
int find(i,j)
int i,j;
{
   return(i*ESTATFIISIZE+j);
}
int tind(i,j)
int i,j;
{
   return(i*ESTATTHETASIZE+j);
}
/*----------------------------------------------------------------------------*/
/*                    Polycrystalline calc. subroutines                       */
/*----------------------------------------------------------------------------*/
getpolysize(poly,reccalc)
struct poly *poly;
struct reccalc *reccalc;
{
   poly->cursize=reccalc->Polysize+gaussianrand()*reccalc->Polysizesigma;
   if (poly->cursize < 0.0) poly->cursize=0.0;
   DEBUGSR("getpolysize:",poly->cursize);
}

exitpolydomain(poly,recoil)
struct poly *poly;
struct recoil *recoil;
{
   poly->prevvel=recoil->vel;
   poly->prevv.x=recoil->vx;
   poly->prevv.y=recoil->vy;
   poly->prevv.z=recoil->vz;
   poly->ndomain++; poly->ndomainstot++;
   poly->boundreached=True;

}


