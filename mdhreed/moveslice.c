#include "common.h"
/*---------------------------------------------------------------------------*/
/*------------------------ SLICE MOVING STUFF -------------------------------*/
/*
  The slice moving algorithms are always tricky to understand, especially
with several layers where the number of atoms may differ. In this 
program the following approach was used:

Remember that there are two different ATOM arrays: 
- The atoms in the simulation cell, whose number _always_ is nat.
  The data of these atoms is stored in structures at,atex and atflag.
- The atoms in the inistate array iniat, whose number may be different for
  each layer. The number in layer nlayer is given by layers->layer[nlayer].natoms

Also remember that movement in several different directions may happen
during one moveslice call ! This has been solved by looping through the entire
subroutine separately for each direction.

Lets label atoms in the simulation cell with atom_s, and atoms in the
inistate array with atom_i

0. The recoil atom (*) has just entered region 4, so:
1. Select an inistate cell nini randomly
2. Select the layer into which the recoil atom has just entered. 
3. Zero the new atom index counter j, and store the recoil atom in a safe place.
4. Go through all atom_s indices i from 0 through nat-1
   5. If the atoms_s position is in regions 2-4 and the atom_s is _not_ 
      the recoil atom, move it backwards by x/4. Give the 
      atom label j. 
      (Since i starts from zero, too, and the recoil atom always has the
      largest index, it is not possible to overwrite the data of a cell in
      regions 2-4. But beware of this in the future !.
   6. If the atom was moved, increase j by 1
7. Go to 4.
8. Go through all atom_i indices from 0 through layers->layer[nlayer].natoms
   9. If the atoms_i inistate position is in region 1, copy it to region 4 in the 
      simulation cell and increase j by 1.
   9.b Call change() to select the atom type randomly
   10. For all inistate cells in layer nlayer, move atoms_i in regions 2-4 backward,
       and atoms_i in region 1 forward.
11. Go to 8
12. Move recoil atom back x/4, and copy it's data into position j. Increase j by 1
13. Set nat = j


 <----- x ----->                 <----- x -----> 
                 
 x/4 x/4 x/4 x/4                 x/4 x/4 x/4 x/4 
 ___ ___ ___ ___                 ___ ___ ___ ___ 
|   |   |   |   |               |   |   |   |   |
|   |   |   |   |               |   |   |   |   |
|   |   | --|*->|  moveslice    |   | --|*->|   |
| 1 | 2 | 3 | 4 |  =========>   | 2 | 3 | 4 | 1'|
|   |   |   |   |               |   |   |   |   |
|   |   |   |   |               |   |   |   |   |
|   |   |   |   |               |   |   |   |   |
|___|___|___|___|               |___|___|___|___|


  Unfortunately, due to the use of structures, x y and z coordinates have to be
handled separately, which results in a lot of extra code. 

The number of atoms in the simulation cell can vary during this process ! 
Therefore, their number nat is recalculated during every call of moveslice. 
If the number of atoms has changed, the number of the recoil atom is also
modified !
*/
/*---------------------------------------------------------------------------*/
moveslice(at,atex,atflag,type,iniat,layers,nlayer,nat,box,
	  recoil,reccalc,recresult,update,firstrecstep,recoilspec,physical)
struct atom *at;
struct atomex *atex;
struct atomflag *atflag;
struct iniat *iniat;
struct type type[MAXTYPE];
struct layers *layers;
int *nlayer,*nat; 
struct box *box;
struct recoil *recoil;
struct reccalc *reccalc;
struct recresult *recresult;
struct physical *physical;
logical *update,firstrecstep,recoilspec;
{
   int i,j,k,l,nini,irec,iii[3];
   struct atom recat;
   struct atomex recatex;
   struct atomflag recatflag;

   struct atom tempat[ORIGMAXATOM];
   struct atomex tempatex[ORIGMAXATOM];
   struct atomflag tempatflag[ORIGMAXATOM];

   int nlayernext,namorphnext;
   int ind,itype;
   int namorph,ll,oldlayer,ii,jj;  /*JARKKO-damage*/
   real E,r,Theta,x,y,z,mx,my,mz;
   point center;

   point move,move1;  /* Holds movement of the recoil atom, -move is movement of the cell */
   logical moved;
   int nmove;
   real xscale,yscale,zscale,rx,ry,rz,tol=0.0001;

   if (! recoil->beeninside) {
     if (recoil->z > box->movelim.z) recoil->beeninside=True; 
        /*ekan layerin ohi, nyt ioni ei valttamatta ole ollut pinnan*/
        /*alla, mutta tarkistus riittaa moveslicen toimintaan. */
   }
   
  
   if(recoil->sz > (box->movelim.z+physical->sput_surf)) recoil->realbeeninside=True; /*aidon pinnan alap.*/

/* Step 0: Check whether movement is needed */
   nmove=0;
   move1.x=move1.y=move1.z=0.0;
   if (recoil->x > box->size.x-box->movelim.x) { move1.x=-box->movelim.x; nmove++; }
   if (recoil->x < box->movelim.x) { move1.x=+box->movelim.x; nmove++; }
   if (recoil->y > box->size.y-box->movelim.y) { move1.y=-box->movelim.y; nmove++; }
   if (recoil->y < box->movelim.y) { move1.y=+box->movelim.y; nmove++; }
   if (recoil->z > box->size.z-box->movelim.z && recoil->beeninside)
      { move1.z=-box->movelim.z; nmove++; }

  
   if (recoil->z < box->movelim.z && 
       ( recoil->realbeeninside || (reccalc->type==2 && recoil->turnedaround)) )
      { move1.z=+box->movelim.z; nmove++; }
   
   if (nmove==0) return;
   if (nmove>3) fprintf(fpdebug,"moveslice Horror Warning: nmove>3 %d\n",nmove);
   *update=True;
   
   oldlayer=*nlayer;    

   for(l=0;l<3;l++) {  /* Select one direction to move in */
      move.x=move.y=move.z=0.0;    
      if (l==0 && move1.x!=0.0) move.x=move1.x;
      else if (l==1 && move1.y!=0.0) move.y=move1.y;
      else if (l==2 && move1.z!=0.0) move.z=move1.z;
      else continue;
      MAXDEBUGSR("moveslice: selected moving direction",l);

      /* Step 1: Select an inistate cell nini randomly */
      nini=unirand(0)*layers->Ninistate;
      if (nini>=layers->Ninistate) nini=layers->Ninistate-1;

      /* Step 2: Select the layer into which the recoil atom has just entered. */
      /*layer+amstate are selected from where the ion has now entered (scaling stays right)..*/
      getlayer(layers,recoil->sz,nlayer,firstrecstep); 
 
      getamorphstate(&namorph,recresult,recoil->sz,box->movelim.z,physical,layers,*nlayer);

      

      /*the cell size can alter during the amorph simulation. box move limits must be scaled too*/
      /*this stage scales the old at-box (size=box.size) to match the new dimensions of the amorphstate*/
      /*(size=[layer][amstate].boxsize).*/
      /*Note that box->movelim must be unitboxsize if this routine is used!!*/

      if(physical->scaling && layers->layer[*nlayer].Namorph>1)
	{
	  xscale=layers->layer[*nlayer].unitboxsize[namorph].x/box->movelim.x;
	  yscale=layers->layer[*nlayer].unitboxsize[namorph].y/box->movelim.y;
	  zscale=layers->layer[*nlayer].unitboxsize[namorph].z/box->movelim.z;
	  
	  if(!(xscale<(1.0+tol) && xscale>(1.0-tol) && 
	       yscale<(1.0+tol) && yscale>(1.0-tol) &&
	       zscale<(1.0+tol) && zscale>(1.0-tol)))
	    {
	      /*if(recoil->sz<40) printf("%lg,%lg,%lg\n",xscale,yscale,zscale);	 */
	    if(xscale>1.2 || yscale>1.2 || zscale>1.2 )
	      printf("LARGE VOLUME CHANGES, MIGHT PRODUCE UNPHYSICAL FEATURES!!\n ");
	    box->movelim.x*=xscale;
	    box->movelim.y*=yscale;
	    box->movelim.z*=zscale;
	    
	    rx=recoil->x;ry=recoil->y;rz=recoil->z; /*recoil coords. before scaling	   */

	    center.x=box->size.x/2;
	    center.y=box->size.y/2;
	    center.z=box->size.z/2;

	    recoil->x=(recoil->x - center.x)*xscale;
	    recoil->y=(recoil->y - center.y)*yscale;
	    recoil->z=(recoil->z - center.z)*zscale;
	  
	    for(i=0; i<*nat; i++)
	      {
	      at[i].x=(at[i].x- center.x)*xscale;
	      at[i].y=(at[i].y- center.y)*yscale;
	      at[i].z=(at[i].z- center.z)*zscale;
	      }
	    box->size.x*=xscale;
	    box->size.y*=yscale;
	    box->size.z*=zscale;
	    move.x*=xscale;
	    move.y*=yscale;
	    move.z*=zscale;
	    recoil->x+=box->size.x/2;
	    recoil->y+=box->size.y/2;
	    recoil->z+=box->size.z/2;
	    rx-=recoil->x;ry-=recoil->y;rz-=recoil->z; /*difference due to scaling*/
	    recoil->ddx+=rx;
	    recoil->ddy+=ry;
	    recoil->ddz+=rz;
	    for(i=0; i<*nat; i++)
	      {
	      at[i].x+=box->size.x/2;/*-rx;  keep the recoil atom steady?*/
	      at[i].y+=box->size.y/2;/*-ry;   */
	      at[i].z+=box->size.z/2;/*-rz;*/
	      }
	    }
	}

      /*..but for the layer to be copied ahead the ion, they are selected with (+movement.z)!!*/
      if(l==2) 
	{
      getlayer(layers,recoil->sz+box->movelim.z,&nlayernext,firstrecstep); 
      getamorphstate(&namorphnext,recresult,recoil->sz+box->movelim.z,box->movelim.z,physical,layers,nlayernext);
	}
      else {nlayernext=*nlayer;namorphnext=namorph;}
      /*printf("namo:%d,nlayer:%d,namonext:%d,nlayernext:%d\n",namorph,*nlayer,namorphnext,nlayernext);*/


      /*move all layers by the move limits.*/
      for(i=0; i<layers->Nlayers; i++)
      for(k=0; k<layers->layer[i].Namorph; k++)
	 {
	 xscale=layers->layer[i].unitboxsize[k].x/box->movelim.x;	 
	 yscale=layers->layer[i].unitboxsize[k].y/box->movelim.y;	 
	 zscale=layers->layer[i].unitboxsize[k].z/box->movelim.z;	 
	 mx=move.x*xscale;
	 my=move.y*yscale;
	 mz=move.z*zscale;
	 /*movelimits for fronts at the same place for each box.*/

        
            if(move.x<0) 
	      {layers->layer[i].front[k].x-=mx;
	      if(layers->layer[i].front[k].x>=layers->layer[i].boxsize[k].x)
	      layers->layer[i].front[k].x-=layers->layer[i].boxsize[k].x;
	      }
	    if(move.x>0) 
	      {layers->layer[i].front[k].x-=mx;
	      if(layers->layer[i].front[k].x<0)
	      layers->layer[i].front[k].x+=layers->layer[i].boxsize[k].x;
	      }
	    if(move.y<0) 
	      {layers->layer[i].front[k].y-=my;
	       if(layers->layer[i].front[k].y>=layers->layer[i].boxsize[k].y)
	       layers->layer[i].front[k].y-=layers->layer[i].boxsize[k].y;
	      }
	    if(move.y>0) 
	      {layers->layer[i].front[k].y-=my;
	      if(layers->layer[i].front[k].y<0)
	      layers->layer[i].front[k].y+=layers->layer[i].boxsize[k].y;
	      }
	    
	    if(i==nlayernext)
	      {
		/*printf("Moving layer %d in z-direction by %lg\n",i,mz);*/
	      if(move.z<0) 
	        {layers->layer[i].front[k].z-=mz;
	        if(layers->layer[i].front[k].z>=layers->layer[i].boxsize[k].z)
	        layers->layer[i].front[k].z-=layers->layer[i].boxsize[k].z;
	        }
	      if(move.z>0) 
	        {layers->layer[i].front[k].z-=mz;
	        if(layers->layer[i].front[k].z<0)
	        layers->layer[i].front[k].z+=layers->layer[i].boxsize[k].z;
	        }	
	    }
	 }/*move layers */

     
	    
	 ii=0;

	 iii[0]=iii[1]=iii[2]=0; /*just for checking which type of atoms are copied..*/


	 /*now we are selecting the coordinates in the next layer*/
       
       xscale=box->movelim.x/layers->layer[nlayernext].unitboxsize[namorphnext].x;
       yscale=box->movelim.y/layers->layer[nlayernext].unitboxsize[namorphnext].y;
       zscale=box->movelim.z/layers->layer[nlayernext].unitboxsize[namorphnext].z;
       if(l!=2 || physical->scaling==0 || (xscale<(1.0+tol) && xscale>(1.0-tol) && 
				   yscale<(1.0+tol) && yscale>(1.0-tol) &&
				   zscale<(1.0+tol) && zscale>(1.0-tol)))
	 {
       for(i=0;i<layers->layer[nlayernext].natoms[namorphnext];i++)
	   {
	   x=iniat[iniind(nlayernext,namorphnext,nini,i)].x;
	   y=iniat[iniind(nlayernext,namorphnext,nini,i)].y;
	   z=iniat[iniind(nlayernext,namorphnext,nini,i)].z;
	   
	   x-=layers->layer[nlayernext].front[namorphnext].x;
	   y-=layers->layer[nlayernext].front[namorphnext].y;
	   z-=layers->layer[nlayernext].front[namorphnext].z;
	   if(x<0) x+=layers->layer[nlayernext].boxsize[namorphnext].x;
	   else if(x>=layers->layer[nlayernext].boxsize[namorphnext].x)
	     x-=layers->layer[nlayernext].boxsize[namorphnext].x;
	   
	   if(y<0) y+=layers->layer[nlayernext].boxsize[namorphnext].y;
	   else if(y>=layers->layer[nlayernext].boxsize[namorphnext].y)
	     y-=layers->layer[nlayernext].boxsize[namorphnext].y;

	   if(z<0) z+=layers->layer[nlayernext].boxsize[namorphnext].z;
	   else if(z>=layers->layer[nlayernext].boxsize[namorphnext].z)
	     z-=layers->layer[nlayernext].boxsize[namorphnext].z;

	   
	   if(x<box->size.x && y<box->size.y && z<box->size.z)	     
	     {
	       copyiniat2(i,ii,iniat,tempat,tempatex,tempatflag,nlayernext,nini,namorphnext,x,y,z);
	       ii++;	        
	       iii[iniat[iniind(nlayernext,namorphnext,nini,i)].type]++;
	     }/*else printf("outside:%f,%f,%f\n",x,y,z);*/

	   }
	 }
       else{ /*scaling needed to attach two layers together*/
	for(i=0;i<layers->layer[nlayernext].natoms[namorphnext];i++)
	   {
	   x=iniat[iniind(nlayernext,namorphnext,nini,i)].x;
	   y=iniat[iniind(nlayernext,namorphnext,nini,i)].y;
	   z=iniat[iniind(nlayernext,namorphnext,nini,i)].z;
	   
	   x*=xscale;
	   y*=yscale;
	   x*=zscale;
   
	   x-=layers->layer[nlayernext].front[namorphnext].x;
	   y-=layers->layer[nlayernext].front[namorphnext].y;
	   z-=layers->layer[nlayernext].front[namorphnext].z;
	   if(x<0) x+=layers->layer[nlayernext].boxsize[namorphnext].x;
	   else if(x>=layers->layer[nlayernext].boxsize[namorphnext].x)
	     x-=layers->layer[nlayernext].boxsize[namorphnext].x;
	   
	   if(y<0) y+=layers->layer[nlayernext].boxsize[namorphnext].y;
	   else if(y>=layers->layer[nlayernext].boxsize[namorphnext].y)
	     y-=layers->layer[nlayernext].boxsize[namorphnext].y;

	   if(z<0) z+=layers->layer[nlayernext].boxsize[namorphnext].z;
	   else if(z>=layers->layer[nlayernext].boxsize[namorphnext].z)
	     z-=layers->layer[nlayernext].boxsize[namorphnext].z;

	   
	   if(x<box->size.x && y<box->size.y && z<box->size.z)	     
	     {
	       copyiniat2(i,ii,iniat,tempat,tempatex,tempatflag,nlayernext,nini,namorphnext,x,y,z);
	       ii++;	        
	       iii[iniat[iniind(nlayernext,namorphnext,nini,i)].type]++;
	     }
	   }
           }
      
      /* Step 3: Zero the new atom index counter j, and store the recoil atom in a safe place. */
      irec=recoil->number;
      copyatom(at+irec,atex+irec,atflag+irec,&recat,&recatex,&recatflag); 
      j=0;
      
      MAXDEBUGSRRR("moveslice: moving",nini,*nlayer,*nat);
      MAXDEBUGSRRR("moveslice am.: moving",nini,*nlayer,namorph); /*JARKKO-damage*/
      MAXDEBUGSRRR("moveslice: moving",move.x,move.y,move.z);
      
      /* Steps 4, 5 and 6 */
      
      for (i=0;i<*nat;i++) {
	 if (i==irec) continue;
	 moved=False;

	 /* Note: if move.x/y/z is 0, no moving is done ! */
	 if (move.x<0 && at[i].x>-move.x) { at[i].x+=move.x; moved=True; }
	 if (move.x>0 && at[i].x<box->size.x-move.x) { at[i].x+=move.x; moved=True; }
	 if (move.y<0 && at[i].y>-move.y) { at[i].y+=move.y; moved=True; }
	 if (move.y>0 && at[i].y<box->size.y-move.y) { at[i].y+=move.y; moved=True; }
	 if (move.z<0 && at[i].z>-move.z) { at[i].z+=move.z; moved=True; }
	 if (move.z>0 && at[i].z<box->size.z-move.z) { at[i].z+=move.z; moved=True; }

	 if (j>i) fprintf(fpdebug,"moveslice Horror Warning: i>j %d %d\n",i,j);
	 if (moved) {
	    MAXDEBUGSRRRRR("mc",i,j,at[i].x,at[i].y,at[i].z);
	    if (i!=j) copyatom(at+i,atex+i,atflag+i,at+j,atex+j,atflag+j); 
	    j++;
	 }
	 else {
	    /* Atom not moved => it will be discarded now */

	    /* Update recoil spectrum array */
	    if (recoilspec) {
	       r=recoil->sz;
	       if (reccalc->Trange == 1) r=recoil->r_proj;
	       if (r > reccalc->recspeczmin && r < reccalc->recspeczmax) {
		  /* Get energy of atom in eV */
		  E=0.5*(at[i].vx*at[i].vx+at[i].vy*at[i].vy+at[i].vz*at[i].vz);

		  Theta=at[i].vz/(sqrt(2.0*E)); 
		  Theta=acos(Theta); 

		  itype=at[i].type;
		  MAXDEBUGSRR("recspec",i,E);
		  recresult->rectot[itype]+=E;
		  recresult->irectot[itype]++;
		  if (E>exp(-0.5/recresult->drecspec)) {
		     ind=(int)(log(E)*recresult->drecspec+0.5);
		     if (ind > RECSPECSIZE) {
		       printf("mdh moveslive ERROR: RECSPECSIZE overflow %d\n",ind);
		       exit(0);
		     }
		     recresult->recspec[itype][ind]+=E;
		     recresult->irecspec[itype][ind]++;
		  }
		 /*(JARKKO) theta-angle*/
		  if(E>reccalc->ETmin && E<=reccalc->ETmax)
		    {
		    ind=(int)(Theta*180/PI);
		    recresult->irecspec_Theta[itype][ind]++;
		    }
		   
		  if (E<1.0) {
		     recresult->rec0[itype]+=E;
		     recresult->irec0[itype]++;
		  }
		 
	       }
	    }
	 }
      } /* Step 7: Loop */
      MAXDEBUGSRR("moveslice: Moved sim. cell atoms",i,j);
      
      /*printf("moved %d out of %d atoms in old box\n",j,*nat);*/

      ll=0;jj=0;
      /* Step 8 */

     
	 for (i=0;i<ii;i++)
	   {
	     moved=False; /*moved means, that atom is copied tempat->at*/
	     if(move.x<0 && tempat[i].x>(box->size.x+move.x)) 
	       {moved=True;}
	     else if(move.x>0 && tempat[i].x<move.x) 
	       {moved=True;}
	     if(move.y<0 && tempat[i].y>(box->size.y+move.y)) 
	       {moved=True;}
	     else if(move.y>0 && tempat[i].y<move.y) 
	       {moved=True;}
	     if(move.z<0 && tempat[i].z>(box->size.z+move.z)) 
	       {moved=True;}
	     else if(move.z>0 && tempat[i].z<move.z) 
	       {moved=True;}
	    if (moved) 
	      if((tempat[i].z+box->movement.z-move.z)>=physical->sput_surf) /*JARKKO-sputtering */
		{
	        copyatom(tempat+i,tempatex+i,tempatflag+i,at+j,atex+j,atflag+j);
	        MAXDEBUGSRRRRR("mc",i,j,at[j].x,at[j].y,at[j].z);
	        change(at,type,j);   /* Step 9.b: select type randomly */
	        j++;jj++;	    
	        } else ll++;
	   }
	 /*if(ll!=0) printf("%d atoms discarded in moveslice 'cause below %g ,boxz= %g ,recz= %g \n",ll,physical->sput_surf,box->movement.z-move.z,recoil->sz);*/
	
	

      /*printf("added %d out of %d atoms to old box (now %d atoms) \n\n",jj,ii,j);*/
      

      MAXDEBUGSRR("moveslice: Moved ini. cell atoms into sim cell",i,j);

      /* Step 12: handle recoil atom */
      recat.x+=move.x; recoil->x=recat.x;
      recat.y+=move.y; recoil->y=recat.y;
      recat.z+=move.z; recoil->z=recat.z;
      copyatom(&recat,&recatex,&recatflag,at+j,atex+j,atflag+j);
      recoil->number=j;
      j++;
      
      for(i=j;i<*nat;i++) zeroatom(at+i,atex+i,atflag+i);
      if (*nat != j) MAXDEBUGSRRR("moveslice: number of atoms changed",*nat,j,*nlayer);


      /* Step 13 */   

      *nat=j;
      box->movement.x-=move.x; box->movement.y-=move.y; box->movement.z-=move.z;

   }  /* End of loop over different directions */


   


}

/*---------------------------------------------------------------------------*/
/*------------------ Small moving subroutines -------------------------------*/

moveineg(coord,lim,box,i,j,l)
real *coord;
real lim,box;
int i,j;
logical *l;
{
   if (*coord < -lim) { *coord+=box+lim; if (i==j) *l=True; }
   else *coord+=lim;
}
moveipos(coord,lim,box,i,j,l)
real *coord;
real lim,box;
int i,j;
logical *l;
{
   if (*coord > box-lim) { *coord-=box-lim; if (i==j) *l=True; }
   else *coord+=lim;
}


copyiniat(i,j,iniat,at,atex,atflag,nlayer,ninistate,namorph)
int i,j;
struct atom *at;
struct atomex *atex;
struct atomflag *atflag;
struct iniat *iniat;
int nlayer,ninistate;
{

   zeroatom(at+j,atex+j,atflag+j);
   /*JARKKO-damage*/
   at[j].x=iniat[iniind(nlayer,namorph,ninistate,i)].x;
   at[j].y=iniat[iniind(nlayer,namorph,ninistate,i)].y;
   at[j].z=iniat[iniind(nlayer,namorph,ninistate,i)].z;
   at[j].vx=iniat[iniind(nlayer,namorph,ninistate,i)].vx;
   at[j].vy=iniat[iniind(nlayer,namorph,ninistate,i)].vy;
   at[j].vz=iniat[iniind(nlayer,namorph,ninistate,i)].vz;
   at[j].type=iniat[iniind(nlayer,namorph,ninistate,i)].type;
   
   atflag[j].moving=True;
}

copyiniat2(i,j,iniat,at,atex,atflag,nlayer,ninistate,namorph,x,y,z)
int i,j;
struct atom *at;
struct atomex *atex;
struct atomflag *atflag;
struct iniat *iniat;
real x,y,z;
int nlayer,ninistate;
{

   zeroatom(at+j,atex+j,atflag+j);
   /*JARKKO-damage*/
   at[j].x=x;/*iniat[iniind(nlayer,namorph,ninistate,i)].x;*/
   at[j].y=y;/*iniat[iniind(nlayer,namorph,ninistate,i)].y;*/
   at[j].z=z;/*iniat[iniind(nlayer,namorph,ninistate,i)].z;*/
   at[j].vx=iniat[iniind(nlayer,namorph,ninistate,i)].vx;
   at[j].vy=iniat[iniind(nlayer,namorph,ninistate,i)].vy;
   at[j].vz=iniat[iniind(nlayer,namorph,ninistate,i)].vz;
   at[j].type=iniat[iniind(nlayer,namorph,ninistate,i)].type;
   
   atflag[j].moving=True;
}

copyatom(at1,atex1,atflag1,at2,atex2,atflag2)
/*
   Copy atom data of one atom from at1 to at2
*/
struct atom *at1,*at2;
struct atomex *atex1,*atex2;
struct atomflag *atflag1,*atflag2;
{
   at2->x  = at1->x;   at2->y  = at1->y;   at2->z  = at1->z;
   at2->vx = at1->vx;  at2->vy = at1->vy;  at2->vz = at1->vz;
   at2->ax = at1->ax;  at2->ay = at1->ay;  at2->az = at1->az;
   at2->a1x= at1->a1x; at2->a1y= at1->a1y; at2->a1z= at1->a1z;
   at2->a2x= at1->a2x; at2->a2y= at1->a2y; at2->a2z= at1->a2z;
   at2->type=at1->type;

   atex2->pot=atex1->pot;
   atex2->force=atex1->force;

   atflag2->moving   = atflag1->moving;
   atflag2->recoil   = atflag1->recoil;
   atflag2->energetic= atflag1->energetic;
   atflag2->interact = atflag1->interact;

}

zeroatom(at,atex,atflag)
/*
   Set atom data to empty values
*/
struct atom *at;
struct atomex *atex;
struct atomflag *atflag;
{
   at->x  = 0.0; at->y  = 0.0; at->z  = 0.0;
   at->vx = 0.0; at->vy = 0.0; at->vz = 0.0;
   at->ax = 0.0; at->ay = 0.0; at->az = 0.0;
   at->a1x= 0.0; at->a1y= 0.0; at->a1z= 0.0;
   at->a2x= 0.0; at->a2y= 0.0; at->a2z= 0.0;
   at->type= 0;

   atex->pot= 0.0;
   atex->force= 0.0;

   atflag->moving   = False;
   atflag->recoil   = False;
   atflag->energetic= False;
   atflag->interact = False;

}











