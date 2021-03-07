#include "common.h"

erodesurface(at,atex,atflag,attemp,atextemp,atflagtemp,recoil,physical,box,nat)
struct atom at[];
struct atomex atex[];
struct atomflag atflag[];
struct atom attemp[];
struct atomex atextemp[];
struct atomflag atflagtemp[];
struct recoil *recoil;
struct physical *physical;
struct box *box;
int *nat;
{
 int j=0,k;
 /*printf("nat: %d ion %d\n",*nat,recoil->number);*/
      for(k=0;k<(*nat-1);k++)
	if(at[k].z>=physical->sput_surf)  /*box.movement.z=0 !!*/
	  {copyatom(&at[k],&atex[k],&atflag[k],&attemp[j],&atextemp[j],&atflagtemp[j]);j++;}
      /*else printf("at %d discarded,z:%g\n",k,at[k].z+box->movement.z);*/
      
      for(k=0;k<j;k++) 
         copyatom(&attemp[k],&atextemp[k],&atflagtemp[k],&at[k],&atex[k],&atflag[k]);

      copyatom(&at[*nat-1],&atex[*nat-1],&atflag[*nat-1],&at[j],&atex[j],&atflag[j]); 
      recoil->number=j;
      j++;

      for(k=j;k<*nat;k++) zeroatom(at+k,atex+k,atflag+k);
     
      printf("%d atoms discarded from mdh.c(ion: %d ) 'cause below %g , recz= %g \n",
	     *nat-j,recoil->number,physical->sput_surf,recoil->sz);
      *nat=j;
}




borders(at,box,nat)
struct atom at[];
struct box *box;
int nat;
{
   int i;

   if (box->periodic) {
      for (i=0;i<nat;i++) {
	 if (at[i].x>box->size.x) at[i].x-=box->size.x;
	 if (at[i].y>box->size.y) at[i].y-=box->size.y;
	 if (at[i].z>box->size.z) at[i].z-=box->size.z;
	 if (at[i].x<0.0) at[i].x+=box->size.x;
	 if (at[i].y<0.0) at[i].y+=box->size.y;
	 if (at[i].z<0.0) at[i].z+=box->size.z;
      }
   }    
   
   if (box->warnoutside) {
      for (i=0;i<nat;i++) {
	 if (at[i].x>box->size.x || at[i].y>box->size.y ||
	     at[i].z>box->size.z || at[i].x<0.0 ||
             at[i].y<0.0 || at[i].z<0.0) 
	    fprintf(fpdebug,"borders.c warning: Atom outside i x y z %d %g %g %g\n",
	       i,at[i].x,at[i].y,at[i].z);
      }
   }    

}

getamorphstate(namorph,recresult,recoilz,boxz,physical,layers,nlayer)
int *namorph;
struct recresult *recresult;
real recoilz,boxz;  /*voi pistaa projektio vaihtoehdot mukaan*/
struct physical *physical;
struct layers *layers;
int nlayer;
{
  int k,i;
  real dope;

  if(physical->Damage) /* && layers->layer[nlayer].Namorph>1) */
    {
      if(recoilz>=0.0) {	
	k=(int) (recoilz/boxz+0.5);
	if(k>=ZBOXES) k=ZBOXES-1;
	else if(k<0) k=0;
	
	dope=recresult->depennucl_box[k]*physical->Dosescaling;
	/*printf("recoilz=%lg,k=%d,dope=%lg\n",recoilz,k,dope);*/
	k=0;  
	for(i=0; i<layers->layer[nlayer].Namorph; i++)
	  {
	    if((layers->layer[nlayer].Dose[i]) > dope) break; 
	    else{k=i;}
	  }
      }else k=0;

  *namorph=k;
  /*  if(recoilz<40) printf("SYV. %g , layer %d ,  Scaled DEPON./SIMKOPPI= %g , VALIT. AMKOPPI %d , eli %g eV/SIMKOPPI\n",recoilz,nlayer,dope,k,layers->layer[nlayer].Dose[k]);*/
    }
  else *namorph=0;

}

/*function that spreads the nuclear deposited energy somehow for damage calculations*/
spread_fdn(r,boxz,recresult,fdn,model)
real r,boxz;
struct recresult *recresult;
real fdn;
int model; 
{
  int k,j,j2;
  real f,V,i,fdng,fraction,dx,it;

  /*update the temp table for future distibution*/
  if(model==1)
    {
    k=(int) (r/boxz+0.5);
    if(k>=ZBOXES) k=ZBOXES-1;
    else if(k<0) k=0;
    recresult->depennucl_boxtemp[k]+=fdn;
    } 
  else if(model==3){ /*gaussian distribution*/
    V=1.0/sqrt(2*PI) * (1/0.5) *5.43; /*constant multiplier, 2=only half of the distribution(area=0.5)*/
     for(k=0;k<ZBOXES; k++)
      {
      i=0;
      f=recresult->depennucl_boxtemp[k];
      if(f>1.0)
	{
	  if(f<100000) i=0.16*pow(f,0.75);  /*after 100kev, lower increase*/
	  else i=1.6*pow(f,0.55);
	             /*some stragling function, fitted in simulated data*/
	/*printf("FD:%g, str:%g \n",f,i);*/
	fraction=1;fdng=0;j=k;
	it=1/(2.0*i*i);
        do{
	  dx=(j-(k+0/5.43)+0.5)*5.43;
	  fraction=V*f/i * exp(-dx*dx*it); 
	                
	  /*if(fraction>100) printf("frac:%g\n",fraction);*/
	  recresult->depennucl_box[j]+=fraction;
	  fdng+=fraction;
	  j++;
	}while(fdng<f && j<ZBOXES && fraction>0.1);  
        if(fdng>f) recresult->depennucl_box[j-1]-=(fdng-f);
	else recresult->depennucl_box[k]+=(f-fdng);
        /*printf("annettu %g , jaettu %g\n",f,fdng);*/
	}
      else recresult->depennucl_box[k]+=f;
      }
  }
  else if(model==4) /*no distribution*/
    {
    k=(int) (r/boxz+0.5);
    if(k>=ZBOXES) k=ZBOXES-1;
    else if(k<0) k=0;
    recresult->depennucl_box[k]+=fdn;
    }
 
}

check_depencalc(physical,elstop,recoil,reccalc,box,poly,recresult,layers,n,Ymean)
struct physical *physical;
struct elstop *elstop;
struct recoil *recoil;
struct reccalc *reccalc;
struct box *box;
struct poly *poly;
struct recresult *recresult;
struct layers *layers;
int n;
real *Ymean;
{
   double dr2,r2,FD2,FDn2,FDe2,Y;
   real ddz=5.43;
   int aa=0;

  if(n==0){  /*if inner loop*/

   if(physical->Damage) 
     {
     FD2=(recoil->eprev-recoil->energy);
     FDe2=elstop->Se*recoil->dr;
     FDn2=FD2-FDe2-elstop->totalfirsov*recoil->dr;

     r2=recoil->sz;  
     if (reccalc->Trange == 1) r2=recoil->r_proj;
     if (reccalc->poly) r2+=poly->range;

     /*energy jakauma boxeittain temp fileen tai suoraan */
     if(physical->spreadFD)
       spread_fdn(r2,box->movelim.z,&(*recresult),FDn2,1);
     else 
       spread_fdn(r2,box->movelim.z,&(*recresult),FDn2,4);
     }

	if(physical->sputtering) 
	   {
	     if(!physical->Damage)
	       {
	       FD2=(recoil->eprev-recoil->energy);
	       FDe2=elstop->Se*recoil->dr;
	       FDn2=FD2-FDe2-elstop->totalfirsov*recoil->dr;

	       r2=recoil->sz;  
	       if (reccalc->Trange == 1) r2=recoil->r_proj;
	       if (reccalc->poly) r2+=poly->range;
	       }
	   /*energy at the surface*/
	   
	   if(r2>=physical->sput_surf && r2<(physical->sput_surf+5.43))
	     {
	     recresult->surfFDn+=FDn2;
	     /*printf("r2= %g Fdn= %g\n",r2,recresult->surfFDn);*/
	     }
	   }
  }
  else{ /*outer loop*/

     if(physical->sputtering==1)
	{ /*sputtering calculation*/
	
	   physical->deltasput_surf=-physical->sput_surf;	  
	   FD2=recresult->surfFDn/5.43;
	   
	   /*templog=False;*/
	   /*getlayer(layers,physical.sput_surf,aa,templog);*/
	   aa=0; /*not exactly right, if the first layer is eaten away..*/
	   printf("ddz: %g layer %d ,deltaFdn: %g \n",ddz,aa,FD2);
	   FDn2=layers->layer[aa].ndensity; /*surface density */
	   if(physical->sputFD==0.0) /*artificial limit needed for spikes*/
	     {	       
	     Y=0.042*FD2/(FDn2*layers->layer[aa].surfbind) - physical->sputdeltaV;
	     *Ymean+=Y;
	     physical->sput_surf+=Y*physical->Dosescaling*5.43/8;	         	      
	     }
	   else 
	     {
	     Y=0.042*physical->sputFD/(FDn2*layers->layer[aa].surfbind) - physical->sputdeltaV;   
	     /*steady surface depon energy*/
	     *Ymean+=Y;
	     physical->sput_surf+=Y*physical->Dosescaling*5.43/8;
	     }  
	   physical->deltasput_surf+=physical->sput_surf; /*erotus edellisia rangeja varten*/
	   /*4.63=E_coh(=surface binding energy for Si[ion-solid-interact. kirja] s.226),
	     7.81ev[handbook of ion impl.tech] ja 2*nnbours/coordinaton*E_coh=1.33*E_coh 
	     [low energy ion irradiation of sol. surf. s91] */	  	    
       }
       else if (physical->sputtering==2)  /*fixded value for sputtering yield*/
	 {
	   physical->deltasput_surf=-physical->sput_surf;
	   *Ymean+=physical->Y - physical->sputdeltaV;
	   physical->sput_surf+=physical->Y*physical->Dosescaling*5.43/8;
	   physical->deltasput_surf+=physical->sput_surf;	  
	 }

 }

}



/*---------------------------------------------------------------------------*/
/*-------- LAYER HANDLING SUBROUTINES ---------------------------------------*/

getlayer(layers,z,nlayer,firstrecstep)
/*
   Subroutine tells in which layer the recoil atom currently is, and returns 
   the value in *nlayer. The various options are designed for maximum
   versability, and must therefore be used with great care: it is quite
   possible to create completely insane combinations of layer selection criteria.

   Note: no checking of multiple layers is done, so the routine
   always selects the last layer that fits the limits. This is
   a bit dangerous, but enables stacking of layers.

   However, layers>1 with limits >< +-1e9 can't be selected.

   If the parameter Probability is 1, and Probmode is 0 (default values) 
   layer will be selected straightly according to minz and maxz, nothing
   special here. If not, things get tricky:
 
   The use of the parameters Probmode and Thickness is:
   Probmode = 0:  Ignore both Probmode and Thickness, select layer randomly at each call.
   Probmode = 1:  Select layer possibly¹ only first time minz is crossed, i.e.
                  only once for each recoil event. If layer is selected, it stays on until
		  the layer limit, unless some other layer with a greater index in the same 
		  region overrides it.
   Probmode = 2:  Select layer possibly¹ every time recoil is between minz and maxz,
                  but once it is selected let it continue atleast for thickness Thickness.
		  When the layer has been thus selected, it stays on until the layer
		  limit, unless some other layer with a greater index in the same region 
		  overrides it.


¹ Possibly here means that layer is selected with probability Probability.

*/
struct layers *layers;
real z;
int *nlayer;
logical firstrecstep;
{ 
   int i,n;
   real unirand();
   static real prevz,thicknessleft=0.0; /* params for Probmode=2 */

   if (firstrecstep) {     /* Zero internal params for new recoil */
      thicknessleft=0.0;
      for (i=0;i<layers->Nlayers;i++) {
	 layers->layer[i].selected=False;
	 layers->layer[i].crossed=False;
      }
   }

   n=-1;

   for (i=0;i<layers->Nlayers;i++) { /* Layers selection loop */
      if (layers->layer[i].minz<=-1e9 && layers->layer[i].maxz>=1e9 && i>0) continue;

      if (z > layers->layer[i].minz && z < layers->layer[i].maxz) { /* Within limits */
	 if (layers->layer[i].Probmode==0  && unirand(0)<=layers->layer[i].Probability) n=i;
	 if (layers->layer[i].Probmode==1) {
	    if (layers->layer[i].crossed && layers->layer[i].selected) { n=i; continue; }
	    else if (layers->layer[i].crossed) continue; 
	    layers->layer[i].crossed=True;
	    layers->layer[i].selected=False;
	    if (unirand(0)<=layers->layer[i].Probability) { n=i; layers->layer[i].selected=True; }
	 }
	 if (layers->layer[i].Probmode==2) {
	    if (thicknessleft>0.0 && thicknessleft<=layers->layer[i].Thickness) {
	       thicknessleft-= z-prevz;
	       MAXDEBUGSRRR("getlayer: ",thicknessleft,z,prevz);
	       n=i;
	    }
	    else if (unirand(0)<layers->layer[i].Probability) { 
	       n=i; 
	       layers->layer[i].selected=True;
	       thicknessleft=layers->layer[i].Thickness;
	    }
	    else layers->layer[i].selected=False;
	 }
      }
   }
   if (n<0) {
      fprintf(fpdebug,"getlayer Warning: No proper layer found at depth %g. Using 0\n",z);
      n=0;
   }
   else MAXDEBUGSRRRRR("getlayer: Layer selected",n,z,layers->layer[n].Probability,layers->layer[n].Probmode,layers->layer[n].Thickness);

   *nlayer=n;
   prevz=z;

}
/*---------------------------------------------------------------------------*/
/*------ Atom changing subroutines ------------------------------------------*/
change(at,type,i)
/*
  If the atom with index i is a changeable atom, select its type
  randomly.
*/
struct atom *at;
struct type type[MAXTYPE];
int i;
{
   real R,probsum;
   int j;
   int oldtype,newtype;

   if(type[at[i].type].prob==1.0 || type[at[i].type].prob==0.0) return;
   oldtype=newtype=at[i].type;

   R=unirand(0);

/* printf("prob i %d type %d prob %g R %g\n",i,at[i].type,type[at[i].type].prob,R); */

   probsum=0.0;
   for (j=0;j<MAXTYPE;j++) {
      if (type[j].prob!=1.0) {
	 probsum+=type[j].prob;
/* printf("prob %d %g %g\n",j,type[j].prob,probsum); */
	 if (R<=probsum) { /* Select this type */
	    newtype=j;
	    break;
	 }
      }
   }
   at[i].type=newtype;
   MAXDEBUGSRRRR("Atom type changed: oldtype newtype R probsum",
		 oldtype,newtype,R,probsum);

   
}

/*---------------------------------------------------------------------------*/
/*------ Miscellaneous small subroutines ------------------------------------*/
real sqdist(x1,y1,z1,x2,y2,z2,box)
/*
   Calculate square of distance between two points, taking border 
   conditions in account if periodic is true.
*/
real x1,y1,z1,x2,y2,z2;
struct box *box;
{
   static real r,r1,r2,r3;

   r1=oabs(x2-x1);
   if (box->periodic && r1+r1>box->size.x) r1=box->size.x-r1;
   r2=oabs(y2-y1);
   if (box->periodic && r2+r2>box->size.y) r2=box->size.y-r2;
   r3=oabs(z2-z1);
   if (box->periodic && r3+r3>box->size.z) r3=box->size.z-r3;
   r=r1*r1+r2*r2+r3*r3;
   return(r);
}

/*---------------------------------------------------------------------------*/
real dist(x1,y1,z1,x2,y2,z2,box)
/*
   Calculate distance between two points, taking border 
   conditions in account if periodic is true.
*/
real x1,y1,z1,x2,y2,z2;
struct box *box;
{
   static real sqr,r,r1,r2,r3;

   r1=oabs((x2-x1));
   if (box->periodic && r1+r1>box->size.x) r1=box->size.x-r1;
   r2=oabs(y2-y1);
   if (box->periodic && r2+r2>box->size.y) r2=box->size.y-r2;
   r3=oabs(z2-z1);
   if (box->periodic && r3+r3>box->size.z) r3=box->size.z-r3;
   r=r1*r1+r2*r2+r3*r3;
   sqr=sqrt(r);
   return(sqr);
}

/*---------------------------------------------------------------------------*/
getmeandisplacement(displ,array,at,nat,box)
real *displ;
point array[];
struct atom *at;
int nat;
struct box *box;
{
   static int i;
   real dist(),r;

   r=0;
   for(i=0;i<nat;i++) {
      r+=dist(at[i].x,at[i].y,at[i].z,array[i].x,array[i].y,array[i].z,box);
   }
   *displ=r/nat;
DEBUGSRRR("getmeandispl",*displ,r,nat);
   return;
}

/*---------------------------------------------------------------------------*/









