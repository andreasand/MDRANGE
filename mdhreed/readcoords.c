#include "common.h"

int borders_for_iniat(coord,lim)
real *coord;
real lim;
{
  int k=0;
   if(*coord < 0) {*coord+=lim;k=1;} 
   else if(*coord>=lim) {*coord-=lim;k=1;}
   if(*coord<0 || *coord>=lim) printf("READCOORDS.C: atom place outside box!!\n");
   return(k);
}

readcoords(iniat,layers,file,type,box,unit,physical)
struct iniat *iniat;
struct layers *layers;
struct file *file;
struct type *type;
struct box *box;
struct unit *unit;
struct physical *physical;
{
   char fname[140],buf[160];
   int i,j,n,nat,t,k,dimflag,rdimflag,y; 
   real a[3],b[3],V,multia,multib,multic; /*JARKKO*/
   int c[3];
   point realboxsize; 
   /*if bscaling, then the coordinate files in layers, that have Namorph>1, have to be*/
   /*origo centered!!!! and the Unit Cell dimensions must be in every file (even with Namorph=0,1)!*/

   printf("NUM OF LAYERS: %d \n",layers->Nlayers);

   for(i=1;i<=layers->Nlayers;i++) {
     
     if(physical->Damage==0 && layers->layer[i-1].Namorph!=1) 
       {
       printf("No damage buildup selected for layer %d, setting namorph=1\n",i);
       layers->layer[i-1].Namorph=1;
       } 
     
     
     for(k=0; k<layers->layer[i-1].Namorph; k++) 
	{
	dimflag=rdimflag=0;
	strcpy(fname,file->coord);
        j=i;
        /*if (layers->layer[i-1].ltype > 0) j=layers->layer[i-1].ltype;*/
	if(physical->Damage==0) /*if no damagebuildup, go with different numbering*/
	  {
	  if (j != 1) 
	    sprintf(fname,"%s.%d",fname,j);  /*coords.in.num_of_layer*/
          }
	else sprintf(fname,"%s.%d.%d",fname,j,k);  
      
      if (debug) fprintf(fpdebug,"Opening coord file %s ltype %d for layer %d am.state %d\n",fname,j,i,k);
      file->fpco=fopen(fname,"r");
      if (file->fpco==NULL) {
	if(debug) fprintf(fpdebug,"readcoords.c: Can't open file %s for layer %d for am state %d\n",fname,i,k);
	if (j==1 && physical->Damage==0) {  /*if damage=0, coords.in naming is ok*/
	    sprintf(fname,"%s.1",fname); 
	    if (debug) fprintf(fpdebug,"Opening coord file %s\n",fname);
	    file->fpco=fopen(fname,"r");
	    if (file->fpco==NULL) 
	       fprintf(fpdebug,"readcoords.c: Still can't open file %s for layer %d for am state %d\n",fname,i,k);

	 }     
      }
      nat=0;
      
      if(k==0) /* only for the first box*/
	{
        layers->layer[i-1].Zmean=0.0;
        layers->layer[i-1].density=0.0;
        layers->layer[i-1].SIdensity=0.0;
	}
      /*boxsize is ALWAYS the real size of the coordinate boxes, even after the scaling*/
      /*set default values first*/
      layers->layer[i-1].boxsize[k]=box->size; 
      /*default values*/
      realboxsize=box->size;
      /*nubmer of boxes for the layers that are not scaled.*/
      c[0]=floor(box->size.x/box->movelim.x+0.5);
      c[1]=floor(box->size.y/box->movelim.y+0.5);
      c[2]=floor(box->size.z/box->movelim.z+0.5);
      multia=multib=multic=1.0;
           
      while(fgets(buf,160,file->fpco)) {
         if (buf[0]=='%' || buf[0]=='!' || buf[0] == '\n') continue;
	 /*read volume as the coordinates are written.*/
	 if(buf[0]=='V') 
	   {
	   n=sscanf(buf+1,"%lf %lf %lf %d %d %d",b+0,b+1,b+2,c+0,c+1,c+2);
	   b[0]*=c[0];b[1]*=c[1];b[2]*=c[2];
	   layers->layer[i-1].boxsize[k].x=b[0];
	   layers->layer[i-1].boxsize[k].y=b[1];
	   layers->layer[i-1].boxsize[k].z=b[2];	   
	  
	   printf("FINISHED READING CELLVOLUME FOR am. state %d: %g,%g,%g\n",k,b[0],b[1],b[2]);
	   dimflag++;
	   
	   if(b[0]<box->Max.x) printf("LAYER SMALLER THAN SIM.BOX!! USING BSCALING??\n");
	   if(b[1]<box->Max.y) printf("LAYER SMALLER THAN SIM.BOX!! USING BSCALING??\n");
	   if(b[2]<box->Max.z) printf("LAYER SMALLER THAN SIM.BOX!! USING BSCALING??\n");
	   continue;
	   }
	 /*if needed, the box can be scaled to match the volume given after "O" */
	 if(buf[0]=='O') 
	   {
	   n=sscanf(buf+1,"%lf %lf %lf",b+0,b+1,b+2);
	   b[0]*=c[0];b[1]*=c[1];b[2]*=c[2];
	   realboxsize.x=b[0];
	   realboxsize.y=b[1];
	   realboxsize.z=b[2];
	   printf("FINISHED READING REAL VOLUME FOR am. state %d: %g,%g,%g\n",k,b[0],b[1],b[2]);
	   rdimflag++;  
	   if(b[0]<box->Max.x) printf("LAYER SMALLER THAN SIM.BOX!! CRITICAL ERROR!!\n");
	   if(b[1]<box->Max.y) printf("LAYER SMALLER THAN SIM.BOX!! CRITICAL ERROR!!\n");
	   if(b[2]<box->Max.z) printf("LAYER SMALLER THAN SIM.BOX!! CRITICAL ERROR!!\n");
	   multia=b[0]/layers->layer[i-1].boxsize[k].x; /*multipliers*/
	   multib=b[1]/layers->layer[i-1].boxsize[k].y;
	   multic=b[2]/layers->layer[i-1].boxsize[k].z;
	   
	   printf("Multipliers read: %lg, %lg, %lg\n",multia,multib,multic);
	   continue;
	   }

	 if(buf[0]=='D') 
	   {
	   n=sscanf(buf+1,"%lf",b+0);
	   layers->layer[i-1].Dose[k]=b[0]*8.0;  /*unit cell*/
	   printf("FINISHED READING DOSE FOR am. state %d: %g(eV/atom)=%g(eV/unitcell)\n",k,b[0],b[0]*8);
	   continue;
	   }
	 
	 n=sscanf(buf,"%lf %lf %lf %d",a+0,a+1,a+2,&t);

         if (n<3) break;
         if (n==3) t=1; 
	 /*JARKKO-damage*/
         iniat[iniind(i-1,k,0,nat)].x=a[0];
         iniat[iniind(i-1,k,0,nat)].y=a[1];
         iniat[iniind(i-1,k,0,nat)].z=a[2];
         iniat[iniind(i-1,k,0,nat)].type=t;
	 file->layeratoms[t][i-1]=1; /*this type (t) of atom is in this layer (i-1)*/
       
          /*if scaling needed in the beginning*/
	 if(physical->bscaling && dimflag && rdimflag)  /*origocentered box!!!*/
	   {
	     /*printf("Doing scaling for atom %d\n",nat);*/
	   	       	       
	   iniat[iniind(i-1,k,0,nat)].x*=multia;
	   iniat[iniind(i-1,k,0,nat)].y*=multib;
	   iniat[iniind(i-1,k,0,nat)].z*=multic;	
	   layers->layer[i-1].boxsize[k]=realboxsize;
	   }
	 
         if (maxdebug) 
	   fprintf(fpdebug,"Coord #%d %g %g %g %d\n",nat,a[0],a[1],a[2],t);
	 if(k==0)
	   {
	   layers->layer[i-1].Zmean+=type[t].Z;	 
	   layers->layer[i-1].density+=1.0;   /* Internal masses are all 1.0 ! */
	   layers->layer[i-1].SIdensity+=unit->mass[t];
	   }
         nat++;
	 if (nat > layers->Natomsmax) {
	    fprintf(fpdebug,"Too many atoms, you should increase layers->Natomsmax from %d\n",layers->Natomsmax);
	    exit(0);
	 }
      } /*while looppi*/
      if(k==0)
	{
        layers->layer[i-1].Zmean/=nat;
V=layers->layer[i-1].boxsize[k].x*layers->layer[i-1].boxsize[k].y*layers->layer[i-1].boxsize[k].z;
        layers->layer[i-1].density/=V;
	layers->layer[i-1].SIdensity/=V*unit->length*unit->length*unit->length;
        layers->layer[i-1].ndensity=nat/(V);
	}
      DEBUG(iniat[0].x);
      if(k==0)
      fprintf(file->fpo,
	      "Read in %d atoms in layer %d amorphization cell %d Zmean %g density %g kg/m^3 i.e. %g atoms/Å^3\n",
	      nat,i,k,layers->layer[i-1].Zmean,layers->layer[i-1].SIdensity,layers->layer[i-1].ndensity);
      else fprintf(file->fpo,"Read in %d atoms in layer %d amorphization cell %d, other data not calculated\n",nat,i,k);

      if (debug) fprintf(file->fpo,"Number density %g at/Å^3\n",layers->layer[i-1].ndensity);
      layers->layer[i-1].natoms[k]=nat;

/*make sure, that all atoms are in a box (0,0,0,sizex,..y,..z)    */
      t=0;
	 for(n=0;n<nat;n++)
	    {
	    y=0;
	    y=borders_for_iniat(&(iniat[iniind(i-1,k,0,n)].x),layers->layer[i-1].boxsize[k].x);
	    y=borders_for_iniat(&(iniat[iniind(i-1,k,0,n)].y),layers->layer[i-1].boxsize[k].y);
	    y=borders_for_iniat(&(iniat[iniind(i-1,k,0,n)].z),layers->layer[i-1].boxsize[k].z);
	    t+=y;
	    }
	 printf("Had to move %d coordinates\n",t);
      /*set borders for layer movement*/
      layers->layer[i-1].front[k].x=0.0;
      layers->layer[i-1].front[k].y=0.0;
      layers->layer[i-1].front[k].z=0.0;
      b[0]=layers->layer[i-1].unitboxsize[k].x=layers->layer[i-1].boxsize[k].x/c[0];
      b[1]=layers->layer[i-1].unitboxsize[k].y=layers->layer[i-1].boxsize[k].y/c[1];
      b[2]=layers->layer[i-1].unitboxsize[k].z=layers->layer[i-1].boxsize[k].z/c[2];
      
      printf("unitboxsize(volume of the coordinates): %lg,%lg,%lg num %d %d %d\n",b[0],b[1],b[2],c[0],c[1],c[2]);
      
	} /*k-loop over am. boxes*/
    
   } /*i-loop over layers*/
   
  
}







