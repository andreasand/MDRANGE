#include "common.h"
/*
  This part is the very heart of the MD simulation code. 
  It is also usually the most CPU-intensive part. Great care should
  be taken in programming it !
*/
/*

  The subroutine works roughly as follows:

  checkdisplacements(updateneeded,updatepos,at,potcrit)
                                      Check displacements from
                    		      previous neighbour list calc,
				      set updateneeded if they're big
				      enough
  getngbrlist(G,updateneeded,updatepos,mode)   
                                      Get neighbour list G and update 
                                      updatepos if updateneeded, 
                                      else skip calculation.

  getaccs(G)   Calculate accelerations for interactions determined in G,
               mark the particles noninteracting if they don't
               experience any interaction.


  The parameter mode determines how the accelerations are to be
  calculated:
  
  mode
  ----
   0         Full calculation of all interactions
   1         Only interactions with recoil atom calculated 
   2         Only interactions between energetic atoms calculated (NOT YET IMPLEMENTED !)

*/
accelerations(at,atex,atflag,pot,potcrit,box,physical,unit,nat,mode,update,file)
struct atom *at;
struct atomex *atex;
struct atomflag *atflag;
struct pot *pot;
struct potcrit *potcrit;
struct box *box;
struct physical *physical;
struct unit *unit;
int nat,mode;
logical *update;
struct file *file;
{
   static int prevmode;
   static logical updateneeded,firsttime=True;
   static point *updatepos;

   static int *G,NG;

   if (!firsttime) updateneeded=False;
   if (firsttime) {
      firsttime=False;
      updateneeded=True;
      MAXNGBR=2*NGBRSTEP;
      updatepos=malloc(MAXATOM*sizeof(point));
      G=(int *) malloc(MAXNGBR*MAXATOM*sizeof(int));
      DEBUGSR("First time in accelerations, allocating ",MAXNGBR*MAXATOM*sizeof(int));
   }

   MAXDEBUGS("Going into checkdisplacements\n");
   if (*update || prevmode != mode || updateneeded) updateneeded=True;
   else checkdisplacements(&updateneeded,updatepos,at,atflag,nat,potcrit,box,mode);
   *update=False;

ngbr:
   MAXDEBUGS("Going into getngbrlist\n");
   getngbrlist(G,&NG,updateneeded,updatepos,mode,at,atflag,nat,potcrit,box);
   if (NG<0) {
      MAXNGBR+=NGBRSTEP;
      DEBUGSRR("Increasing neighbour list size to ",MAXNGBR*MAXATOM*sizeof(int),NG);
      if (debug) fprintf(fpdebug,"%p\n",G);
      G=realloc(G,MAXNGBR*MAXATOM*sizeof(int));
      if (debug) fprintf(fpdebug,"%p\n",G);
      NG=0;
      goto ngbr;
   }
   MAXDEBUG(NG);
   
   MAXDEBUGS("Going into getaccs\n");
   getaccs(G,NG,at,atex,atflag,pot,potcrit,box,physical,unit,nat,mode,file);

   prevmode=mode;
}

/*---------------------------------------------------------------------------*/
checkdisplacements(updateneeded,updatepos,at,atflag,nat,potcrit,box,mode)
logical *updateneeded;
point updatepos[];
struct atom *at;
struct atomflag *atflag;
int nat;
struct potcrit *potcrit;
struct box *box;
int mode;
{
   static real r1,r2,r3,r12,r22,r32,r2lim;
   static real halfx,halfy,halfz,sizex,sizey,sizez;
   static int i;

   r2lim=(potcrit->rm-potcrit->r0)*(potcrit->rm-potcrit->r0);

   halfx=box->size.x/2.0; sizex=box->size.x;
   halfy=box->size.y/2.0; sizey=box->size.y;
   halfz=box->size.z/2.0; sizez=box->size.z;

   for(i=0;i<nat;i++) {
      if (mode==1 && !atflag[i].recoil) continue;
      r1=oabs(updatepos[i].x-at[i].x);
      if (box->periodic && r1>halfx) r1=sizex-r1;
      r12=r1*r1;
      if (r12 > r2lim) goto end;
      r2=oabs(updatepos[i].y-at[i].y);
      if (box->periodic && r2>halfy) r2=sizey-r2;
      r22=r2*r2;
      if (r12+r22 > r2lim) goto end; 
      r3=oabs(updatepos[i].z-at[i].z);
      if (box->periodic && r3>halfz) r3=sizez-r3;
      r32=r3*r3;
      if (r12+r22+r32 < r2lim) continue; 

end:  *updateneeded=True; 
       break; 
   }
}
/*---------------------------------------------------------------------------*/
getngbrlist(G,NG,updateneeded,updatepos,mode,at,atflag,nat,potcrit,box)
/*
   The neighbour list G contains atom indices. Because 0 is a valid atom
index, the value -1 is used to denote that a new atom should be selected,
and the value -2 denotes the end of the list.

*/
int *G,*NG;
logical updateneeded;
point updatepos[];
int mode;
struct atom *at;
struct atomflag *atflag;
int nat;
struct potcrit *potcrit;
struct box *box;
{
   real r2lim;
   int i,i0,j;
   real r,r1,r2,r3,r12,r22,r32;
   real halfx,halfy,halfz,sizex,sizey,sizez;
   int fa=0;
   if (!updateneeded) return; 
   MAXDEBUGS("Doing neighbour list update\n");
   r2lim=potcrit->rm*potcrit->rm;
         
   *NG=0;
   halfx=box->size.x/2.0; sizex=box->size.x;
   halfy=box->size.y/2.0; sizey=box->size.y;
   halfz=box->size.z/2.0; sizez=box->size.z;
   
   for(i=0;i<nat;i++) {
      G[*NG]=-1; (*NG)++;
      updatepos[i].x=at[i].x;
      updatepos[i].y=at[i].y;
      updatepos[i].z=at[i].z;

      if (mode==1 && !atflag[i].recoil) continue;
      if (mode==2 && !atflag[i].recoil && !atflag[i].energetic) continue;

      i0=i+1;
      if (mode==1 || mode==2) i0=0;
      for (j=i0;j<nat;j++) {
	 if (i==j) continue;
 	 r1=oabs(at[j].x-at[i].x);
	 if (box->periodic && r1>halfx) r1=sizex-r1;
	 r12=r1*r1;
	 if (r12 > r2lim) continue;
	 r2=oabs(at[j].y-at[i].y);
	 if (box->periodic && r2>halfy) r2=sizey-r2;
	 r22=r2*r2;
	 if (r12+r22 > r2lim) continue; 
	 r3=oabs(at[j].z-at[i].z);
	 if (box->periodic && r3>halfz) r3=sizez-r3;
	 r32=r3*r3;
	 if (r12+r22+r32 > r2lim) continue; 

	 G[*NG]=j; (*NG)++;
	 if(i==nat-1) fa++;
MAXDEBUGSRRR("nl",i,j,r12+r22+r32);
	 if (*NG>=MAXNGBR*MAXATOM) { *NG=-1; return; }
       }
   } 
   G[*NG]= -2;  /* End of list */
   if (maxdebug && *NG<=nat) {
      fprintf(fpdebug,"accelerations.c Warning: no interactions\n");
      printco(at,nat,box,0.0,0);
   }
   /*printf("%d naapuria at Z:%g,radius:%g\n",fa,at[nat-1].z+box->movement.z,potcrit->rm);*/
   /*if((at[nat-1].z+box->movement.z)<200 && fa!=0 ) printf("z: %g nbrs: %d\n",at[nat-1].z+box->movement.z,fa);*/
   MAXDEBUGSR("acc: ",*NG); 
}
/*---------------------------------------------------------------------------*/
getaccs(G,NG,at,atex,atflag,pot,potcrit,box,physical,unit,nat,mode,file)
int G[],NG;
struct atom *at;
struct atomex *atex;
struct atomflag *atflag;
struct pot *pot;
struct potcrit *potcrit;
struct box *box;
struct physical *physical;
struct unit *unit;
int nat,mode;
struct file *file;
{
   int i,j,c;
   real r,r2lim,V,dV;
   real dVr,dV1,dV2,dV3;
   real r1,r2,r3,r12,r22,r32;
   real halfx,halfy,halfz,sizex,sizey,sizez;
   logical recoil;
   real ir; /*JARKKO*/

   for (i=0;i<nat;i++) {
      atflag[i].interact=False;
      at[i].ax=at[i].ay=at[i].az=0.0;
      atex[i].pot=0.0;
      atex[i].force=0.0;
   }
   potcrit->fatoms=0; /*firsov neighbours*/

   physical->epot=0.0;
   r2lim=(potcrit->r0)*(potcrit->r0);
   halfx=box->size.x/2.0; sizex=box->size.x;
   halfy=box->size.y/2.0; sizey=box->size.y;
   halfz=box->size.z/2.0; sizez=box->size.z;

   c= -1; i= -1;     /* i=-1 because first element of G is always -1 */
   while(True) {     /* Used to be while (c<NG) */
      c++;
      j=G[c];
      if (j==-2) break;
      if (j==-1) { i++; continue; }  /* Previous atom done, start next atom neighbour loop */
      r1=at[j].x-at[i].x;
      if (box->periodic) { if (r1>halfx) r1-=sizex; else if (r1<-halfx) r1+=sizex; }
      r12=r1*r1;
      if (r12 > r2lim) continue;
      r2=at[j].y-at[i].y;
      if (box->periodic) { if (r2>halfy) r2-=sizey; else if (r2<-halfy) r2+=sizey; }
      r22=r2*r2;
      if (r12+r22 > r2lim) continue;
      r3=at[j].z-at[i].z;
      if (box->periodic) { if (r3>halfz) r3-=sizez; else if (r3<-halfz) r3+=sizez; }
      r32=r3*r3;
      if (r12+r22+r32 > r2lim) continue;
      
      r=sqrtf(r12+r22+r32);
      atflag[i].interact=atflag[j].interact=True;
      
      recoil=atflag[i].recoil||atflag[j].recoil;
      ir=1/r;
      potential(ir,r,&V,&dV,at[i].type,at[j].type,recoil,pot,potcrit,file);
      
      dVr=dV*ir; dV1=dVr*r1; dV2=dVr*r2; dV3=dVr*r3;
      
      /*firsov neighbours*/
      if(r<potcrit->rfirsov) 
	if(potcrit->fatoms<MAXN){
           potcrit->frmin[potcrit->fatoms]=r;
           potcrit->fin[potcrit->fatoms]=j;
           potcrit->fatoms++;
           }

      
      at[i].ax+=dV1;
      at[i].ay+=dV2;
      at[i].az+=dV3;
      
      at[j].ax-=dV1;
      at[j].ay-=dV2;
      at[j].az-=dV3;

      atex[i].force+=dV;   /* This doesn't take in account that the force is a vector ! */
      atex[j].force+=dV;
      physical->epot+=V;
      atex[i].pot+=V;
      atex[j].pot+=V;
   }

}

