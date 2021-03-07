#include "common.h"

#define NPOTARRAY 1000

potential(ir,r,V,dV,typei,typej,recoil,pot,potcrit,file)
real ir,r,*V,*dV;
int typei,typej;
logical recoil;
struct pot *pot;
struct potcrit *potcrit;
struct file *file;
{

   static logical firsttime=True;
   static logical firstwarn=False;
   static int i,j;
   static real dr,x,y,y2;
   static double rr,xd,yd;
   static double Vrd,dVrd;
   static real Vr,dVr,Va,dVa,Vj,dVj; /* repulsive, attractive and joined pot. */

   int t1,t2,jmax;
   static real rmax=10.0,step;
   static real repV[MAXTYPE][MAXTYPE][NPOTARRAY],
     repdV[MAXTYPE][MAXTYPE][NPOTARRAY];


/*
   Morse parameters. Parameter 0 means the values are read in from the parameter
file. Other values of pot->attr.morsetype mean:
                      ^^^^^^^^^^^^^^^^^^^
 1   2   3   4   5   6 
 --  --  --  --  --  --
 Ta  Nb  Fe  Ni  W   Cr

Note that Ta and Nb are not based on physical data, other than the
nearest-neighbour distance ! Also, 8.12 1994 D was modified so that the 
potential yields the rms thermal displacement calculated from the Debye
model.
*/
/*                                           Ta    Nb      Fe      Ni      W       Cr    */
   static real morsealpha[MAXMORSE] = { 0.0, 0.9,  1.0, 1.3885, 1.4419, 1.4116, 1.5721 };
   static real morser1[MAXMORSE] =    { 0.0,2.87, 2.86, 2.845,  2.780,  3.032,  2.754 };
   static real morseD[MAXMORSE] =     { 0.0, 2.2,  1.75, 0.4174, 0.4205, 0.9906, 0.4414 };


   step=rmax/NPOTARRAY;

   if (firsttime) {
      firsttime=False;
      if (debug) fprintf(fpdebug,"First time in potential\n");

      DEBUGS("Copying morse parameters\n");
      for (i=1;i<MAXMORSE;i++) {     /* Note the 1, array index 0 is readin ! */
	 pot->attr.morsealpha[i]=morsealpha[i];
	 pot->attr.morser1[i]=morser1[i];
	 pot->attr.morseD[i]=morseD[i];
      }

      DEBUGS("Storing screening function in array\n");

      jmax=rmax/step;
      for (j=0;j<jmax;j++) for (t1=0;t1<MAXTYPE;t1++) for (t2=0;t2<MAXTYPE;t2++) {
	 rr=(double) j*step;
         Vrd=dVrd=0.0;
	 xd=rr/pot->rep.au[t1][t2];
	 for (i=0;i<MAXPOTEXP;i++) {
            if (pot->rep.ai[i]==0.0) break;
	    yd=pot->rep.a[t1][t2][i]*exp(xd*pot->rep.b[t1][t2][i]);
	    Vrd+=yd;
	    dVrd+=pot->rep.bi[i]*yd;
	 }
	 repV[t1][t2][j] =Vrd;
	 repdV[t1][t2][j]=dVrd;
	 
	 
	 if (j==0) MAXDEBUGSRRRR("reppot j=0:",t1,t2,Vrd,dVrd);
      }
      DEBUGSRR("Storing finished",jmax,rmax);


   }  /* End of if firsttime */


   if (recoil) {

      if (pot->rep.type[typei][typej]==1) {

	 /* Repulsive potential of ZBL type (but exponential terms may differ) */
         Vrd=dVrd=0.0;
	 xd=r/pot->rep.au[typei][typej];
	 for (i=0;i<MAXPOTEXP;i++) {
            if (pot->rep.ai[i]==0.0) break;
	    yd=pot->rep.a[typei][typej][i]*exp(xd*pot->rep.b[typei][typej][i]);
	    Vrd+=yd;
	    dVrd+=pot->rep.bi[i]*yd;
	 }

	 *V=pot->rep.scale*Vrd*pot->rep.fact[typei][typej]*ir;
	 *dV=pot->rep.scale*(pot->rep.fact[typei][typej]* dVrd/pot->rep.au[typei][typej] - *V)*ir;	 
	 MAXDEBUGSRRR("pot",r,*V,*dV);
      }

      /* Repulsive potential from spline interpolation */

      if (pot->rep.type[typei][typej]==2) {
	 splinereppot(ir,r,V,dV,typei,typej,pot->rep.fact[typei][typej],file);
      }

      if (pot->rep.type[typei][typej]==3) {

	 /* Repulsive potential of ZBL type from array */
	 if (r>=rmax-step) { 
	    Vrd=dVrd=0.0;
	 }
	 else {
	    j=(int) floor(r/step);
	    dr=r-step*j;
	    Vrd =((step-dr)*repV[typei][typej][j]+dr*repV[typei][typej][j+1])/step;
	    dVrd =((step-dr)*repdV[typei][typej][j]+dr*repdV[typei][typej][j+1])/step;
	 }

	 *V=pot->rep.scale*Vrd*pot->rep.fact[typei][typej]*ir;
	 *dV=pot->rep.scale*(pot->rep.fact[typei][typej]* dVrd/pot->rep.au[typei][typej] - *V)*ir;	 


	 DEBUGSRRR("pot",r,*V,*dV);
      }


      MAXDEBUGSRRR("pot",r,*V,*dV);
      return;

   }


/* Attractive potential */

   Va=dVa=0.0;
   if (pot->attr.type == 1) {   /* Mazzone attractive */
      if (r < pot->attr.mazr2-pot->attr.mazd) {
         y=expf(-pot->attr.mazalpha*(r-pot->attr.mazr1));
         y2=y*y;
         Va=pot->attr.mazD*(y2-y-y);
         dVa=pot->attr.mazD*pot->attr.mazalpha*(-y2+y);
         dVa+=dVa;
      }
      else if (r > pot->attr.mazr2-pot->attr.mazd && 
               r < pot->attr.mazr2+pot->attr.mazd) {
         y=r-pot->attr.mazr2;
         Va=pot->attr.mazK*y*y;
         dVa=pot->attr.mazK*y;
         dVa+=dVa;
      }
   }

   if (pot->attr.type == 2) {   /* Morse */
      i=pot->attr.morsetype;
      y=expf(-pot->attr.morsealpha[i]*(r-pot->attr.morser1[i]));
      y2=y*y;
      Va=pot->attr.morseD[i]*(y2-y-y);
      dVa=pot->attr.morseD[i]*pot->attr.morsealpha[i]*(-y2+y);
      dVa+=dVa;
   }

   if (pot->attr.type < 1 || pot->attr.type > 2) 
           fprintf(fpdebug,"potential.c: invalid potential type !\n");

   *V=Va;
   *dV=dVa;
  /* printf("po %g %g %g\n",r,*V,*dV); */
}





/**************************************************************************/

/**************************************************************************/

/**************************************************************************/



splinereppot(ir,r,V,dV,ti,tj,screenfact,file)
real ir,r;         /* input r */
real *V,*dV;    /* Output result */
int ti,tj;      /* Atom types i and j */
real screenfact; /* Coulombic scale factor Z1 Z2/4 pi eps_0 */ 
struct file *file;
{
/*
   Spline repulsive potential subroutine, taken from ~/md/potential/tersoff.c

   - first time each pair is requested, read in potential and derivatives from reppot.in.ti.tj

   - always: return V and dV for all atoms i

   - interpolation of V and dV done with cubic splines. See Numerical
     Recipes, ch. 3.3

   - Due to the linear extrapolation to 0 (?), the spline interpolation tends to
     give a somewhat bad derivative between the 2nd and 3rd points. The error
     is not large, and we usually are not that high in energy anyway.

   - Screening version: input scaled by screenfact/r for spline interpolation,
     result scaled back again. Thus the interpolated factor is really the
     screening function, which is less strongly repulsive and thus better
     suited for spline interpolation

     screenfact should always contain Z1 Z1/4 pi eps_0 in the correct
     units, i.e. 14.399758 * Z1 Z2 in units of eV Å.

     The screening factor calc. was tested for Kr->Pt up to an energy
     of 10 MeV an seemed to work perfectly - no nasty spline oscillations.
*/

   const real rmax=9.0;             /* max r, if r bigger than rmax set V to zero */

   static real x[REPPOTMAX][MAXTYPE][MAXTYPE];                   /* x array */ 
   static real VR[REPPOTMAX][MAXTYPE][MAXTYPE];                  /* Repulsive potential */ 

   static real VR2[REPPOTMAX][MAXTYPE][MAXTYPE];    /* Estimation of second derivative of VR for spline int. */

   /* Spline interpolation temporary arrays and help variables */

   static real u[REPPOTMAX];       /* Spline help array */

   static int N;           /* Number of read in potential points */
   static int klo=0,khi=1;  

   real sig,p,qn,un;
   int k;

   real h,a,b,rlo,rhi;
   real screenhelp,dVspline;

   /* Other variables */

   static int firsttime=True;
   static int firstpair[MAXTYPE][MAXTYPE];

   char buf[120],*s;
   int i,j,l,n;
   double rin,Vin,rinprev;

   real x1,x2,x3,y1,y2,y3;

   real r0,deltar;

   
   real temp,atemp,btemp,htemp;


   MAXDEBUGSRRR("splinereppot()",r,ti,tj);

   /* On firsttime read in potential data */
   
   if (firsttime) {
      firsttime=False;
      for (i=0;i<MAXTYPE;i++) {
	 for (j=0;j<MAXTYPE;j++) {
	    firstpair[i][j]=True;
	 }
      }
   }

   if (firstpair[ti][tj]) {

      firstpair[ti][tj]=False;
      sprintf(buf,"reppot.%d.%d.in",ti,tj);

      file->fpre=fopen(buf,"r");
      if (file->fpre==NULL) {
	 fprintf(file->fpdb,"splinereppot() ERROR: You requested read in reppot for atom type pair %d %d\n",
                         ti,tj);
	 fprintf(file->fpdb,"But I can't open file %s\n\n",buf);
	 exit(0);
      }

      fprintf(file->fpo,"Reading in reppot for atom type pair %d %d, r %g\n",ti,tj,r);

      i=0; l=0;
      rinprev=0.0;
      while (True) {
	 s=fgets(buf,120,file->fpre);
	 if (s==NULL) break;

	 if (buf[0]=='%' || buf[0]=='#' || buf[0]=='!' || buf[0] == '\n') continue;

	 l++;
	 n=sscanf(buf,"%lg %lg",&rin,&Vin);
	 if (n>=2) {
	    i++;
	    if (i>=REPPOTMAX) fprintf(file->fpdb,"reppot() ERROR: too many reppot.%d.%d points %d\n\n",ti,tj,l);
	    x[i][ti][tj]=rin;
	    VR[i][ti][tj]=Vin/screenfact*rin;
	    if (rin <= rinprev) fprintf(file->fpdb,"reppot() Warning: x data in wrong order on line %d\n\n",l);
	    rinprev=rin;
	    
	 }
	 else {
	    fprintf(file->fpdb,"reppot() Warning: too few arguments on line %d\n\n",l);
	 }
      }
      N=i;
      fprintf(file->fpdb,"reppot() Read in %d valid-seeming repulsive potential points\n",N);
      fclose(file->fpre);

      /* Fill rest of arrays with r:s and 0 */
      for(j=N+1;j<REPPOTMAX;j++) {
	 x[j][ti][tj]=x[N][ti][tj]+1.0*(j-N)*(x[N][ti][tj]-x[N-1][ti][tj]);
	 VR[j][ti][tj]=0.0;

	 VR2[j][ti][tj]=0.0;
	 u[j]=0.0;
      }

      /* 
	 Obtain 0 point estimate.
	 Since we deal with the screening factor, we know the exact values !
      */
      x[0][ti][tj]=0.0;
      VR[0][ti][tj]=1.0;

      /*
	 Initialize spline stuff for spline interpolation of VR 

	 Input needed: vector VR, size REPPOTMAX
	 Note that array bounds are 0 and N,
	 not 1 and N ! 

	 Output: VR2, estimate of second derivative of VR (corresponds
	 to Y2 in Numerical Recipes)

      */

      VR2[0][ti][tj]=u[0]=0.0;
      
      for (i=1;i<=N-1;i++) {
	 sig=(x[i][ti][tj]-x[i-1][ti][tj])/(x[i+1][ti][tj]-x[i-1][ti][tj]);
	 p=sig*VR2[i-1][ti][tj]+2.0;
	 VR2[i][ti][tj]=(sig-1.0)/p;
	 u[i]=(6.0*((VR[i+1][ti][tj]-VR[i][ti][tj])/
		    (x[i+1][ti][tj]-x[i][ti][tj]) - 
		    (VR[i][ti][tj]-VR[i-1][ti][tj]) 
		    /(x[i][ti][tj]-x[i-1][ti][tj]))/
	       (x[i+1][ti][tj]-x[i-1][ti][tj])-sig*u[i-1])/p;
      }
      
      qn=un=0.0;

      VR2[N][ti][tj]=(un-qn*u[N-1])/(qn*VR2[N-1][ti][tj]+1.0);
      for (k=N-1;k>=0;k--) {
	 VR2[k][ti][tj]=VR2[k][ti][tj]*VR2[k+1][ti][tj]+u[k];
      }
	  

      fprintf(file->fpo,"Spline initialization of screening functions done for atom type pair %d %d\n",ti,tj);
   
   } /* End of if (firstpair[ti][tj]) */


   *V=0.0;
   *dV=0.0;
   
   if (r>=rmax || r>=x[N][ti][tj]) return;
	    
   /* Obtain array limits klo and khi */

   if ((khi-klo == 1) && x[khi][ti][tj] > r && x[klo][ti][tj] < r) ;
   else { /* bisect your way to the correct point */
      klo=0; 
      khi=N;
      while (khi-klo > 1) {
	 k=(khi+klo)/2;
	 if (x[k][ti][tj] > r) khi=k;
	 else klo=k;
      }
   }

   rlo=x[klo][ti][tj];
   rhi=x[khi][ti][tj];

   /* Spline interpolate V */

   h=rhi-rlo;
   temp=1/h;
   a=(rhi-r)*temp;
   b=(r-rlo)*temp;
   
   atemp=a*a;
   btemp=b*b;
   htemp=h*h/6.0;
   
   *V=a*VR[klo][ti][tj]+b*VR[khi][ti][tj]+((atemp*a-a)*VR2[klo][ti][tj]+(btemp*b-b)*VR2[khi][ti][tj])*htemp;

   /* Compute spline derivative, eq. 3.3.5 in Numerical recipes */

   dVspline=((VR[khi][ti][tj]-VR[klo][ti][tj])-
     ((3*atemp-1)*VR2[klo][ti][tj]+(3*btemp-1)*VR2[khi][ti][tj])*htemp)*temp;

   /* Correct for screening; */

   screenhelp=screenfact*ir;
   *dV=screenhelp*(r*dVspline- *V)*ir;
   *V  *= screenhelp;

/*
   *V/=2.0;
   *dV/=2.0;
*/

   /* printf("p %g %g %g\n",r,*V,*dV); */

}


