#include "common.h"

readelstop(file,elstop,type,layers)
struct file *file;
struct elstop *elstop;
struct type *type;
struct layers *layers;
{
   int i,n,j,counts,k,l,jakaja;
   char buf[170],fname[140];
   real a,b,alfa;
   char *p;

   elstop->penr=0;   

  
   for(j=0; j<layers->Nlayers; j++)  /*layer loop*/
     {
       printf("type of elstop for layer %d is %d\n",j,layers->layer[j].elstop);
       if(layers->layer[j].fulldensity)
	 {
	 printf("Allocating %d bytes for charge density table in layer %d\n",
		MAXRHO*MAXRHO*MAXRHO*sizeof(real),j);
	 
	 elstop->rho=(real *)malloc(MAXRHO*MAXRHO*MAXRHO*sizeof(real));
	 
	 counts=0;
	 printf("Reading charge density file for full box..");   
         file->fprho=fopen(file->crho,"r");
	 if(file->fprho==NULL) printf("Can't read %s\n",file->crho);

	 
	 for(i=0;i<MAXRHO;i++)
	 for(l=0;l<MAXRHO;l++)
	 for(k=0;k<MAXRHO;k++){
	    fscanf(file->fprho,"%lg",&alfa);		 
	    elstop->rho[i+(l*MAXRHO)+(k*MAXRHO*MAXRHO)]=alfa;
	    counts++;
	    }
	 fclose(file->fprho);
	 printf("DONE, %d points\n",counts);
	 }

       /*elstop type loop*/
       if(layers->layer[j].elstop==0)  /*max number of zbl-layers=ORIGMAXLAYER*/
	 {
	 strcpy(fname,file->elstop);
	 if(layers->Nlayers>1)  /*if multiple layers, elstop.in.[1,2,3..], else elstop.in*/
           sprintf(fname,"%s.%d",fname,j+1);
	 file->fpel=fopen(fname,"r");
	 if (file->fpel == NULL) fprintf(fpdebug,"Can't read elstop(ZBL) file %s\n",fname);
	 if (file->fpel == NULL && layers->Nlayers<2) /*only one layer, but no elstop.in, trying elstop.in.0 */
	   {
	   sprintf(fname,"%s.%d",fname,j);
	   fprintf(fpdebug,"Trying to open elstop(ZBL) file %s\n",fname);
	   file->fpel=fopen(fname,"r");
	   if (file->fpel == NULL) fprintf(fpdebug,"Can't read elstop(ZBL) file %s\n",fname);
	   }

	 i=0;
	 while(fgets(buf,160,file->fpel) != NULL) {
	   if (buf[0]=='%' || buf[0]=='!' || buf[0] == '\n') continue;
	   n=sscanf(buf,"%lf %lf",&a,&b);
	   if (n<2) break;
	   elstop->v[j][i]=a; elstop->Setable[j][i]=b;
	   i++;
	   if(i >= MAXELSTOP) 
	     {
	     fprintf(fpdebug,"Too many points in elstop array, trying to manage with less, increase MAXELSTOP\n");
	     fflush(fpdebug);
	     break;
	     }
	   }
	 elstop->n[j]=i;
	 fprintf(file->fpo,"Read in %d elstop(ZBL) points from %s for layer %d , max vel %g max Se %g\n",
		 elstop->n[j],fname,j,elstop->v[j][i-1],elstop->Setable[j][i-1]);
	 }

       /* puska or BK, same el.density file*/
       if(layers->layer[j].elstop==1 || layers->layer[j].elstop==2)
	 {
	   elstop->penr++;
         if(layers->layer[j].fulldensity==0) /*if not full density file provided*/
	   {
	     /*this checks if the ion type is also occuring in the layers!!*/
	     l=1;
	     for(k=0; k<ORIGMAXLAYER; k++)
	       if(file->layeratoms[0][k]!=0) l=0;

	     layers->layer[j].steadydens=0;
	     jakaja=0; /*counts the number of types in layer, needed in steadydens.*/

	     for(k=l;k<MAXTYPE;k++) /*from 1, if type 0 is only the ion, from 0 if not !*/
	       if(file->layeratoms[k][j]!=0) /*if this type (k) of atom is in this layer (j)*/
	       {
	       strcpy(fname,file->cdens);	       
               sprintf(fname,"%s.%d",fname,k);
 
	       file->fpcdens=fopen(fname,"r");
	       if(file->fpcdens == NULL)
		 {
	         fprintf(fpdebug,"Can't read charge density file %s\n",fname);
		 strcpy(fname,file->cdens);
		 fprintf(fpdebug,"Trying to open file %s\n");
		 file->fpcdens=fopen(fname,"r");
		 if(file->fpcdens == NULL)
		   {
		   fprintf(fpdebug,"Can't read charge density file %s\n",fname);
		   fprintf(fpdebug,"No charge density for atom type %d in the layer %d, hopefully its OK?\n",k,j);
		   goto jumpty;
		   }
		 else 
		   fprintf(fpdebug,"Using %s for atom type %d in layer %d, hopefully its OK?\n",fname,k,j);
		 }
	       i=0;
	       while(fgets(buf,160,file->fpcdens) != NULL) {
	         if (buf[0]=='%' || buf[0]=='!' || buf[0] == '\n') continue;
	         n=sscanf(buf,"%lg %lg",&a,&b);
	         if (n<2) break;
	         elstop->cdist[k][i]=a; elstop->celdens[k][i]=b;
	         i++;
	         if(i >= MAXSTOPP) 
	           {
	           fprintf(fpdebug,"Too many points in cdens array, trying to manage with less, increase MAXSTOPP\n");
	           fflush(fpdebug);
	           break;
	           }
	         }    
	       elstop->nc[k]=i;
	       layers->layer[j].steadydens+=elstop->celdens[k][i-1];
	       jakaja++;
	       fprintf(file->fpo,"Read in %d cdensity points from %s for layer %d atomtype %d, last density=%g\n",
		       elstop->nc[k],fname,j,k,elstop->celdens[k][i-1]);
	       jumpty:;
	       }
	   if(jakaja>0) layers->layer[j].steadydens/=jakaja;
	   else layers->layer[j].steadydens=0;
	   fprintf(file->fpo,"Steady el.dens in layer %d is %lg\n",j,layers->layer[j].steadydens);
	   }
   
   /*puska, stopping.in is Se versus r_s*/
   if(layers->layer[j].elstop==1) 
     {
     file->fpstopp=fopen(file->stopp,"r");
     if (file->fpstopp == NULL) fprintf(fpdebug,"Can't read stopping vs. rs file %s\n",file->stopp);
 
     i=0;
     while(fgets(buf,160,file->fpstopp) != NULL) {
      if (buf[0]=='%' || buf[0]=='!' || buf[0] == '\n') continue;
      n=sscanf(buf,"%lf %lf",&a,&b);
      if (n<2) break;
      elstop->strs[i]=a; elstop->stopp[i]=b;
        if(i==1) elstop->slim=(elstop->strs[1]-elstop->strs[0]);
        if(i>1) {
        if(100*fabs(elstop->slim-elstop->strs[i]+elstop->strs[i-1])>elstop->slim) elstop->slim=0;
        }
      i++;
      if (i >= MAXSTOPP) {
         fprintf(fpdebug,"Too many points in stopping vs. rs array, trying to manage with less, increase MAXSTOPP\n");
         fflush(fpdebug);
         break;
      }
   }
   elstop->ns=i;
   fprintf(file->fpo,"Read in %d el.stopping points from %s for layer %d , last stopping=%g\n",
	   elstop->ns,file->stopp,j,elstop->stopp[i-1]);
     }   
   else if(layers->layer[j].elstop==2) 
     {
       /*BK, elstop->stopp[] is now G*S_lin versus r_s, which we precalculate now*/
       b=6.0/MAXSTOPP;  /*maximum r_s is 6.0 for G (jussi's thesis)*/
       alfa=pow((4.0/(9.0*PI)),0.3333333);
       for(i=0; i<MAXSTOPP; i++)
	 {
	 a=(i+1)*b;
	 elstop->strs[i]=a;
	 elstop->stopp[i]=2.0/(3.0*PI)*( log(1.0+PI/(alfa*a))-1.0/(1.0+alfa*a/PI) );
	 elstop->stopp[i]*=1.0+0.717*a-0.125*a*a-0.0124*a*a*a+0.00212*a*a*a*a;
	 elstop->stopp[i]*=2.34951e-05; /* unit conversion atomic units-> ev/ong */
	 if(i==1) elstop->slim=(elstop->strs[1]-elstop->strs[0]);
         if(i>1) 
           if(100*fabs(elstop->slim-elstop->strs[i]+elstop->strs[i-1])>elstop->slim) elstop->slim=0;
	 /*printf("%lg %lg\n",a,elstop->stopp[i]);*/
	 }     
     elstop->ns=MAXSTOPP;
     fprintf(file->fpo,"Precalculated %d el. stopping points for layer %d , last value=%g\n",
	   elstop->ns,j,elstop->stopp[i-1]);
     }

   if(elstop->slim!=0) elstop->iperslim=1.0/elstop->slim;  /*needed in the future, just a speedup JARKKO*/
   else elstop->iperslim=1.0;

   /*DEBUGSR("ZEFF", elstop->zlim);*/
   printf("Distance between stopping density points=%lg\n", elstop->slim);
  
       } /*puska,BK if*/
     } /*layer loop*/
 
}

real getelstop(elstop,unit,recoil,v,nlayer)
struct elstop *elstop;
struct unit *unit;
struct recoil *recoil;
int nlayer;
real v;
{
   int i,ihalf,n;
   int iup,idown;
   real vreal,Se;
   real scale;
   real omega,LSSOmega();

   vreal=v*unit->vel[recoil->attype];
   n=elstop->n[nlayer];
   /* 
      Get two nearest elements with binary search 
   */
   iup=n-1; idown=0;
   while (True) {
      ihalf=idown+(iup-idown)/2;
      if (ihalf==idown) break;
      if (elstop->v[nlayer][ihalf] < vreal) idown=ihalf;
      else iup=ihalf;
   }

   /* Calculate scaling factor using the formula scale+dscale*E+sqscale*E*E */
   scale=elstop->scale+recoil->energy*(elstop->dscale+elstop->sqscale*recoil->energy);

   /* Get elstop with a simple linear interpolation */
   Se=elstop->Setable[nlayer][idown]*(elstop->v[nlayer][iup]-vreal)/(elstop->v[nlayer][iup]-elstop->v[nlayer][idown]);
   Se+=elstop->Setable[nlayer][iup]*(vreal-elstop->v[nlayer][idown])/(elstop->v[nlayer][iup]-elstop->v[nlayer][idown]);
   Se*=scale;

   return(Se);

}

subelstop(at,time,recoil,elstop,unit,type,layer,potcrit,nlayer,file)
struct atom *at;
struct time *time;
struct recoil *recoil;
struct elstop *elstop;
struct unit *unit;
struct type *type;
struct layer layer;  /* Current layer ! */
struct potcrit *potcrit;
int nlayer ; /*nuber of current layer*/
struct file *file;
{
   int i,ihalf,n,k;
   static int iup=0,idown=1;
   real v,vreal,deltav,Se,Semod;
   real scale;
   real omega,LSSOmega();
   int jup=0,jdown=1,jjhalf;   
   int jjup=0,jjdown=1,ty;
   real temp,charge,onelrad,kf,ktf;
  
   real z1,z2,alpha,ACR,BCR,SCR,fvx,fvy,fvz,firtotal=0,avx,avy,avz;  /*firsov*/
   real firtotalx=0,firtotaly=0,firtotalz=0;
   int ain,atyp;
   real gamma,q,vr,vf; /*BK*/
   real ui,vi,ti,Se2,plasmon_f;
   double elx,ely,elz,rex,rey,rez;
   int xup,yup,zup,xdown,ydown,zdown;
   int axd,axu,ayd,ayu; /*(JARKKO)*/
   real *apurho; /*(JARKKO)*/
   real ti_ui,ti_vi,vi_ui,ti_ui_vi; /*(JARKKO)*/
   
   
   elstop->Se=Se=charge=0.0;

   i=recoil->number;
   v=recoil->vel;
   vreal=v*unit->vel[recoil->attype];


   if(layer.elstop==0) /*ZBL model for this layer?*/
     {      
     if(debug) printf("Using zbl elstop at depth %g layer %d \n",recoil->sz,nlayer);
     n=elstop->n[nlayer];
     
     if (vreal > elstop->v[nlayer][n-1]) {
       printf("\nElstop array maximum exceeded in ZBL-file : vion %g vmaxelstop %g\n\n",
	    vreal,elstop->v[nlayer][n-1]);
      exit(0);
      }

     if (!(vreal > elstop->v[nlayer][idown] && vreal < elstop->v[nlayer][iup])) {
      /* If the previous borders aren't valid, get two nearest 
         elements with binary search 
      */
      iup=n-1; idown=0;
      while (True) {
	 ihalf=idown+(iup-idown)/2;
	 if (ihalf==idown) break;
	 if (elstop->v[nlayer][ihalf] < vreal) idown=ihalf;
	 else iup=ihalf;
         }
      }
     /* Calculate scaling factor using the formula scale+dscale*E+sqscale*E*E */
     scale=elstop->scale+recoil->energy*(elstop->dscale+elstop->sqscale*recoil->energy);

     /* Get elstop with a simple linear interpolation */
     Se=( elstop->Setable[nlayer][idown]*(elstop->v[nlayer][iup]-vreal)
	 +elstop->Setable[nlayer][iup]*(vreal-elstop->v[nlayer][idown])
	 )/(elstop->v[nlayer][iup]-elstop->v[nlayer][idown]);

     Se*=scale;

/*  
   Include straggling of electronic stopping if mode is 1
   See A. Luukkainen and M. Hautala, Radiat. Eff. 59 (1982)
*/
     if (elstop->mode==1 && recoil->dr>0.0) {
      /* printf("elstop %g",Se); fflush(NULL);*/
      omega=LSSOmega(type[recoil->attype].Z,layer.Zmean,recoil->vel,
		     layer.ndensity,recoil->dr,unit,recoil->attype);
      Semod=gaussianrand()*omega/recoil->dr;
      MAXDEBUGSRRR("LSS mod",Se,Semod,omega);
      Se+=Semod; 
      }
   
   }   
   else if(layer.elstop==1)   /*puska model*/
     {
     if(debug) printf("using Puska elstop at depth %g ",recoil->sz);
  
     if(layer.fulldensity==0)
       {
       if(potcrit->fatoms==0) 
	 {
	 if(debug) printf("\nPuska model, no atoms inside rfirsov at depth %g, using default charge density ",recoil->sz); 
	 goto jump2;
	 }
	      
       for(k=0;k<potcrit->fatoms;k++)  /*Si-atoms inside Firsov (default=1.47) radius*/
        {
	  
        ty=at[potcrit->fin[k]].type;
	  if(!(potcrit->frmin[k] > elstop->cdist[ty][jjdown] && potcrit->frmin[k] <elstop->cdist[ty][jjup]))
	    {
	      /* If the previous borders(0-1) aren't valid, get two nearest
		 elements with binary search
	      */
	    jjup=elstop->nc[ty]-1; jjdown=0;
	    while (True) {
	      jjhalf=jjdown+(jjup-jjdown)/2;
	      if(jjhalf==jjdown) break;
	      if(elstop->cdist[ty][jjhalf] < potcrit->frmin[k]) jjdown=jjhalf;
	        else jjup=jjhalf;
	      }
	    }
	/* Get charge density with a simple linear interpolation */
        /* Varaustiheys */
	  /*temp=1.0/(elstop->cdist[jjup]-elstop->cdist[jjdown]);*/
	  charge+=(elstop->celdens[ty][jjdown]*(elstop->cdist[ty][jjup]-potcrit->frmin[k])+elstop->celdens[ty][jjup]*(potcrit->frmin[k]-elstop->cdist[ty][jjdown]))/(elstop->cdist[ty][jjup]-elstop->cdist[ty][jjdown]);
	  /*SP+=elstop->celdens[jjup]*(potcrit->frmin[k]-elstop->cdist[jjdown])*temp;*/
	  MAXDEBUGSRR("elstop.c", Se, potcrit->frmin[k]);	 	 
	} /*fatoms looppi*/

       /*       printf("SP:%g ",SP);*/

jump2:;
       if(charge>0) 
	 {charge*=a_b*a_b*a_b;}  
         /*unit conversion: 1 e/Ongs^3 -> atomic units=1 e/(Bohrradius=0.529 Ong)^3 */
       else 
	 {charge=layer.steadydens*a_b*a_b*a_b;}           
           

       

       onelrad=pow(3.0/(4*PI*charge),0.33333333);  /*one-electron radius from el. charge density      */
       }
     else{
    
        rex=recoil->x;
	rey=recoil->y;
	rez=recoil->z;

	if(rex<=elstop->ux)  /*For silicon ux=5.43 ong */
	  elx=rex;		
	else
	  elx=(rex-(floor(rex/elstop->ux)*elstop->ux));
		
        if(rey<=elstop->uy) 
          ely=rey;  	
        else     
          ely=(rey-(floor(rey/elstop->uy)*elstop->uy));
                 
	/*rez=recoil->sz-potcrit->oxide; */
        if(rez<=elstop->uz) 
          elz=rez;          
        else     
          elz=(rez-(floor(rez/elstop->uz)*elstop->uz));        

	if(elz<0 || elx<0 || ely<0) {
	  printf("cdens variable(%lg,%lg,%lg) LESS THAN ZERO z%lg sz%lg\n",elx,ely,elz,recoil->z, recoil->sz);	
	  exit;
	}
	
	  xup=(int)ceil(elx/elstop->dx);
	  if(xup>0) xdown=xup-1;
	    else    xdown=0;
	  if(xup==128)	xup=0;
	  yup=(int)ceil(ely/elstop->dy);
	  if(yup>0) ydown=yup-1;
	    else    ydown=0;
	  if(yup==128) yup=0;         
	  zup=(int)ceil(elz/elstop->dz);
	  if(zup>0) zdown=zup-1;
	    else    zdown=0;
	  if(zup==128) zup=0;
	  
	  if(zup<0 || zdown>127 || xup<0 || xdown>127 || yup<0 || ydown>127 )
	    {
	    printf("OUT OF BOUNDS ERROR %d %d %d\n",xdown,ydown,zdown); 
	    printf("%lg %lg %lg\n",elx,ely,elz);
	    exit;
 	    }
	  ti=(elx-(xdown*elstop->dx))/elstop->dx;
	  ui=(ely-(ydown*elstop->dy))/elstop->dy;
      	  vi=(elz-(zdown*elstop->dz))/elstop->dz;
	  
	  
	  ayd=ydown<<7;
	  ayu=yup<<7;  
	  axd=xdown<<14;
	  axu=xup<<14;
	  apurho=(real *)(elstop->rho);
	  ti_ui=ti*ui;
	  ti_ui_vi=ti_ui*vi;
	  ti_vi=ti*vi;
	  vi_ui=vi*ui;
	  temp=(1-vi-ti+ti_vi-ui+vi_ui+ti_ui-ti_ui_vi)*apurho[axd+ayd+zdown]
	         +(ti-ti_ui-ti_vi+ti_ui_vi)*apurho[axu+ayd+zdown];
temp+=(ti_vi-ti_ui_vi)*apurho[axu+ayd+zup]+(vi-ti_vi-vi_ui+ti_ui_vi)*apurho[axd+ayd+zup];
temp+=(ui-ti_ui-vi_ui+ti_ui_vi)*apurho[axd+ayu+zdown]+(ti_ui-ti_ui_vi)*apurho[axu+ayu+zdown];
temp+=ti_ui_vi*apurho[axu+ayu+zup]+(vi_ui-ti_ui_vi)*apurho[axd+ayu+zup];

 onelrad=exp(temp);  /*note that rs was in reality log(rs) in atomic units!! */
	  }/*charge if */


     if(elstop->correction)
       if(onelrad>3.0) onelrad*=pow(5.0,-0.33333); /*efective r_s for low el. densities in channels (n<0.06 a.u.)*/
      
     /*printf("one el.rad:%g\n",onelrad);*/

     
       if(onelrad>0.1) {	
         jup=(int)ceil(onelrad/elstop->slim);
         jdown=(int)floor(onelrad/elstop->slim);
         jdown--;
         jup--;
         if(jup>=elstop->ns){
           printf("TOO LOW ELECTRON DENSITY number:%d, using MAX=%d, onelrad:%g\n",jup,elstop->ns,onelrad);
           jup=jdown=elstop->ns-1;
           }
	 else if(/*jdown<0 ||*/ jup < 0){
          printf("TOO HIGH ELECTRON DENSITY number:%d, using 0\n",jdown);
          jdown=jup=0;
          }
	if(jup!=jdown){
	  Se=(elstop->stopp[jdown]*(elstop->strs[jup]-onelrad) + elstop->stopp[jup]*(onelrad-elstop->strs[jdown]))*elstop->iperslim;  /*iperslim= 1/slim*/
	  /*Se+=elstop->stopp[jup]*(onelrad-elstop->strs[jdown])*elstop->iperslim;*/
	  }
	else Se=elstop->stopp[jdown];	
        }
        else { /*Linear response theory*/
	  Se=type[recoil->attype].Z*type[recoil->attype].Z*0.21220659*(log(1+6.0292139/onelrad)-(1/(1+0.16585911*onelrad)));
	  Se*=2.34951e-05; /*atomic units-> [eV] [Ong] [m/s]*/
	}       

       Se*=vreal;
      
     }
   else if(layer.elstop==2)  /*BK-stopping,ala. Cai-Beardmore style*/
     {
     if(layer.fulldensity==0)
       { 

      if(potcrit->fatoms==0) 
	 {
	 if(debug) 
	   printf("BK model, no atoms inside rfirsov at depth %g, using default charge density ",recoil->sz); 
	 goto jump2bk;
	 }

      
       if(debug) printf("using BK elstop at depth %g ",recoil->sz);
       for(k=0;k<potcrit->fatoms;k++)  /*Si-atoms inside Firsov (default=1.47) radius*/
          {
	  ty=at[potcrit->fin[k]].type;  
	  if(!(potcrit->frmin[k] > elstop->cdist[ty][jjdown] && potcrit->frmin[k] <elstop->cdist[ty][jjup]))
	    {	    
	    jjup=elstop->nc[ty]-1; jjdown=0;
	    while (True) {
	      jjhalf=jjdown+(jjup-jjdown)/2;
	      if(jjhalf==jjdown) break;
	      if(elstop->cdist[ty][jjhalf] < potcrit->frmin[k]) jjdown=jjhalf;
	        else jjup=jjhalf;
	      }
	    }
	
	  charge+=(elstop->celdens[ty][jjdown]*(elstop->cdist[ty][jjup]-potcrit->frmin[k])+elstop->celdens[ty][jjup]*(potcrit->frmin[k]-elstop->cdist[ty][jjdown]))/(elstop->cdist[ty][jjup]-elstop->cdist[ty][jjdown]);	 
	  MAXDEBUGSRR("elstop.c", Se, potcrit->frmin[k]);	 	 
	} /*fatoms loop*/

      

jump2bk:;
       
       if(charge>0) 
	 {charge*=a_b*a_b*a_b;/*printf("Sein: %g\n",charge);*/}  
         /*unit conversion: 1 eV/ong^3 -> 1 e/(Bohrradius=0.529 Ong)^3 =0.14803589*/
       else 
	 {charge=layer.steadydens*a_b*a_b*a_b;}           
              

       onelrad=pow(3.0/(4*PI*charge),0.33333333);  /*one-electron radius from el. charge density*/
       }
      else{ 
     /* Let's interpolate the charge density */
        rex=recoil->x;
	rey=recoil->y;
	rez=recoil->z;

	if(rex<=elstop->ux)  /*For silicon ux=5.43 ong */
	  elx=rex;		
	else
	  elx=(rex-(floor(rex/elstop->ux)*elstop->ux));
		
        if(rey<=elstop->uy) 
          ely=rey;  	
        else     
          ely=(rey-(floor(rey/elstop->uy)*elstop->uy));
                 
	/*rez=recoil->sz-potcrit->oxide; */
        if(rez<=elstop->uz) 
          elz=rez;          
        else     
          elz=(rez-(floor(rez/elstop->uz)*elstop->uz));        

	if(elz<0 || elx<0 || ely<0) {
	  printf("cdens variable(%lg,%lg,%lg) LESS THAN ZERO z%lg sz%lg\n",elx,ely,elz,recoil->z, recoil->sz);	
	  exit;
	}
	
	  xup=(int)ceil(elx/elstop->dx);
	  if(xup>0) xdown=xup-1;
	    else    xdown=0;
	  if(xup==128)	xup=0;
	  yup=(int)ceil(ely/elstop->dy);
	  if(yup>0) ydown=yup-1;
	    else    ydown=0;
	  if(yup==128) yup=0;         
	  zup=(int)ceil(elz/elstop->dz);
	  if(zup>0) zdown=zup-1;
	    else    zdown=0;
	  if(zup==128) zup=0;
	  
	  if(zup<0 || zdown>127 || xup<0 || xdown>127 || yup<0 || ydown>127 )
	    {
	    printf("OUT OF BOUNDS ERROR %d %d %d\n",xdown,ydown,zdown); 
	    printf("%lg %lg %lg\n",elx,ely,elz);
	    exit;
 	    }
	  ti=(elx-(xdown*elstop->dx))/elstop->dx;
	  ui=(ely-(ydown*elstop->dy))/elstop->dy;
      	  vi=(elz-(zdown*elstop->dz))/elstop->dz;
	  
	  
	  
	  ayd=ydown<<7;
	  ayu=yup<<7;  
	  axd=xdown<<14;
	  axu=xup<<14;
	  apurho=(real *)(elstop->rho);
	  ti_ui=ti*ui;
	  ti_ui_vi=ti_ui*vi;
	  ti_vi=ti*vi;
	  vi_ui=vi*ui;
	  onelrad=(1-vi-ti+ti_vi-ui+vi_ui+ti_ui-ti_ui_vi)*apurho[axd+ayd+zdown]
	         +(ti-ti_ui-ti_vi+ti_ui_vi)*apurho[axu+ayd+zdown];
onelrad+=(ti_vi-ti_ui_vi)*apurho[axu+ayd+zup]+(vi-ti_vi-vi_ui+ti_ui_vi)*apurho[axd+ayd+zup];
onelrad+=(ui-ti_ui-vi_ui+ti_ui_vi)*apurho[axd+ayu+zdown]+(ti_ui-ti_ui_vi)*apurho[axu+ayu+zdown];
onelrad+=ti_ui_vi*apurho[axu+ayu+zup]+(vi_ui-ti_ui_vi)*apurho[axd+ayu+zup];

 onelrad=exp(onelrad);  /*note that charge was in reality log(rs)!! */
	  }/*charge if */

    

     if(elstop->correction)
       if(onelrad>3.0) onelrad*=pow(5.0,-0.33333); /*efective r_s for low el. densities in channels (n<0.06 a.u.)*/


     /*calculating gamma: Z=Z_eff=gamma*Z1*/
     if(type[recoil->attype].Z>1)
       {           
        vf= 1.9191583/onelrad;   
	vr=vreal/2187640;
	if(vr<vf)
	  q=0.75*vf*(1+(2*vr*vr/(3*vf*vf))-(vr*vr*vr*vr/(15*vf*vf*vf*vf)));
        else
	  q=vf*(1+(vf*vf/(5*vr*vr)));

	/*q is now the relative velocity*/

	q*=1/pow(1.0*type[recoil->attype].Z,0.66667); /*q is now reduced relative vel.*/
	q=1-exp(-0.95*(q-0.07)); /*q is now charge fraction*/
        /*q=0.89*q-0.3*q*q;	 Klein et al., APL 57 , neglible difference(jussi's thesis)*/
        gamma=2*0.24005*pow((1-q),0.66667)/(pow(1.0*type[recoil->attype].Z,0.33333)*(1-((1-q)/7))); /*BK:n screening length lambda*/
	kf=onelrad*pow((9.0*PI)/4,0.33333);
	ktf=sqrt(4*kf/PI);
	gamma=gamma/(1.0-(2.0/3.0)*ktf*ktf*gamma*gamma); /*Kaneko, improvement Prev. A41, 4889 (1990)*/
        gamma=q+0.5*(1-q)*log(1+(16*gamma*gamma/(onelrad*onelrad)));
        gamma*=type[recoil->attype].Z;    
        gamma*=gamma;
        }
     else gamma=1;

      jup=(int)ceil( onelrad/elstop->slim);
      jdown=(int)floor( onelrad/elstop->slim);
      jup--;
      jdown--;
      if(jup>=elstop->ns){
	printf("TOO LOW DENSITY %d\n",jup);
	jup=jdown=elstop->ns-1;
	}
      if(jdown<0){
        printf("TOO HIGH DENSITY %d\n",jdown);
        jdown=jup=0;
        }    
 
      /*Now elstop->stopp must be the correction factor G(r_s)*/
      if(jup!=jdown){
	  Se=elstop->stopp[jdown]*(elstop->strs[jup]-onelrad)/elstop->slim;
	  Se+=elstop->stopp[jup]*(onelrad-elstop->strs[jdown])/elstop->slim;}
      else
	  Se=elstop->stopp[jdown];
      
     Se*=vreal;		
     /* Calculate scaling factor using the formula scale+dscale*E+sqscale*E*E */
     scale=elstop->scale+recoil->energy*(elstop->dscale+elstop->sqscale*recoil->energy);
     Se*=scale; /* We scale the stopping of proton, default value: scale=1 */
     Se*=gamma; /* and multiply it by the square of the effective charge to get the stopping of a heavy ion (gamma=1 for hydrogen)*/
     
     }

   if(file->printdens)
     {
     charge=3.0/(4.0*PI*onelrad*onelrad*onelrad);    
     k =(int) (charge/(10.0/DEPENSIZE)+0.5); 
     /*printf("%d %i\n",k,k);*/
     if(k>=DEPENSIZE) k=DEPENSIZE-1;
     else if(k<0) k=0;
     elstop->dens[k]+=1;
     elstop->Sespc[k]+=Se;
     }

   

    if(debug) printf("S_e is %g \n",Se);
   
   deltav=time->dt[recoil->attype]*Se;
   temp=1-deltav/v;
   at[i].vx*=temp; 
   at[i].vy*=temp;
   at[i].vz*=temp;

   /*Firsov inelastic electronic stopping due to the change of electrons between nuclei*/
   if(elstop->Firsov){
   if(type[recoil->attype].Z>1){
     /* We don't use the Firsov model for protons */
     avx=(at[i].vx*unit->vel[recoil->attype]);
     avy=(at[i].vy*unit->vel[recoil->attype]);
     avz=(at[i].vz*unit->vel[recoil->attype]);
     for(k=0;k<potcrit->fatoms;k++) /* A loop over all colse neighbours */{
       /* For more info on the Firsov parametrization used, see Elteckov et al. in Atomic 
	  Collision Phenomena in Solids, p. 657 */  
       ain=potcrit->fin[k];
       atyp=at[ain].type;
       if(type[recoil->attype].Z > type[at[ain].type].Z){
	 z1 = type[recoil->attype].Z;
	 z2 = type[at[ain].type].Z;
       }
       else {
	 z2 = type[recoil->attype].Z;
	 z1 = type[at[ain].type].Z;
       }  
       alpha=1/(1+pow((double)z2/z1, 0.166666667));
       ACR=(1+0.8*alpha*pow((double)z1,0.33333333)*potcrit->frmin[k]/0.46850229);
       ACR=pow(ACR,4.0);
       BCR=(1+0.8*(1-alpha)*pow((double)z2,0.33333333)*potcrit->frmin[k]/0.46850229);
       BCR=pow(BCR,4.0);
       SCR=(z1*z1/ACR + z2*z2/BCR)*5.240844e-06;  /*5.240844e-06=gfac*/
       fvx=avx-(at[potcrit->fin[k]].vx*unit->vel[at[potcrit->fin[k]].type]);
       fvy=avy-(at[potcrit->fin[k]].vy*unit->vel[at[potcrit->fin[k]].type]);
       fvz=avz-(at[potcrit->fin[k]].vz*unit->vel[at[potcrit->fin[k]].type]);
       firtotalx+=SCR*fvx;   
       firtotaly+=SCR*fvy;
       firtotalz+=SCR*fvz;
       
       
       /*at[potcrit->fin[k]].vx+=SCR*fvx*time->dt[at[potcrit->fin[k]].type];*/
       /*at[potcrit->fin[k]].vy+=SCR*fvy*time->dt[at[potcrit->fin[k]].type];*/
       /*at[potcrit->fin[k]].vz+=SCR*fvz*time->dt[at[potcrit->fin[k]].type];*/
       
       /*firtotal+=sqrt(firtotalx*firtotalx+firtotaly*firtotaly+firtotalz*firtotalz);*/
     }
    
     firtotal=sqrt(firtotalx*firtotalx+firtotaly*firtotaly+firtotalz*firtotalz);
     /* printf("tot firsov %g\n",firtotal);*/
     elstop->totalfirsov=firtotal;
     /*firtotal*=time->dt[recoil->attype];  Force -> change in velocity, dv = (F/m) *dt  
     at[i].vx*=1-(firtotal/v);
     at[i].vy*=1-(firtotal/v);   
     at[i].vz*=1-(firtotal/v);*/ 
     at[i].vx-=firtotalx*time->dt[recoil->attype];
     at[i].vy-=firtotaly*time->dt[recoil->attype];
     at[i].vz-=firtotalz*time->dt[recoil->attype];   
   }
   }else elstop->totalfirsov=0;
   elstop->Se=Se;

   MAXDEBUGSRRR("elstop.c",elstop->v[nlayer][idown],elstop->v[nlayer][iup],Se);

}   

real LSSOmega(Z1,Z2,v,N,deltar,unit,t)
int Z1;
real Z2,v,N,deltar;
struct unit *unit;  
int t;       /* Type of recoil atom */
/*
   Calculate the Omega of the LSS model (Lindhard, Scharff, Schiott, 
   Mat. Fys. Medd. Dan. Vid. Selsk. 33, 14 (1963))
*/
{
   real x,omegaB,L;
   real Bohrvel=2187700.0;   /* SI units */
   real e=1.602189e-19,eps0=8.854188e-12;    
   real vSI;

   omegaB=(Z1*e/(4*PI*eps0)*sqrt(4*PI*Z2*N*deltar/(unit->length*unit->length)));

   vSI=v*unit->vel[t]; 
   x=vSI*vSI/(Bohrvel*Bohrvel*Z2);

   L=1.0;
   if (x < 2.3) {
      L=0.5*(1.36*sqrt(x)-0.016*sqrt(x*x*x));
   }

   MAXDEBUGSRRRR("LSSOmega",v,x,L,omegaB*sqrt(L));

   return(omegaB*sqrt(L));

}


