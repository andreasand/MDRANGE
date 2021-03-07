#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(int argv,char* argc[]) {

int n, npts,nx,ny,nz,m,n1s,n3s,n2s,n2p,n3p;
FILE *fp,*fp2;
FILE *gp;
char ch;
double r, phi1, phi2,phi2p, phi3,phi3p;
float xx,x0,dx,yy,y0,dy,zz,z0,dz;
float y[4][1000];
double s1[8]={0.968,0.03033,0.02248,-0.00617,0.00326,-0.00143,0.00081,-0.0016};
double s2[8]={-0.25755,-0.00446,0.11153,0.40339,0.55032,0.03381,-0.00815,0.00126};
double s3[8]={0.06595,0.00185,-0.03461,-0.10378,-0.19229,-0.06561,0.59732,0.55390};
double s3p[8]={-0.11535,-0.00189,-0.00473,-0.07552,0.01041,0.46075,0.57665,0.06274};
double s2p[8]={0.54290,0.00234,0.04228,0.32155,0.22474,0.00732,-0.00105,0.00041};
double ksiexp[8]={14.01420,16.39320,10.8795,7.72709,5.165,2.97451,2.14316,1.31306};
double ksiexp3p[8]={7.14360,16.25720,10.79720,6.89724,4.66598,2.32046,1.33470,0.79318};
double ksi[8];
double ksi3p[8];
double ksi0[8];
double ksi2p[8], ksi2s[8];
double nfac[8];
double nfac3p[8];
double Y00=0.28209479;
double ksisum, totsum;
double kappa3[2]={0.9847779,0.82582169};
double pi=3.14159265358979323846;
double x2, y2, z2;
double V=160.0/(27.0*sqrt(3.0)*pi);  //constant in (4b)
double a_b=0.529177343575629789; //Ongs
double ax,ay,az,r2,r3,r4;
double hartree=27.2116; //eV

double unitcellx=5.43,unitcelly=5.43,unitcellz=5.43;

//parametrization, n=4 in Deutsch, PRB45,646 (1992)
 double kappa[2]={0.9949,0.9382}; //+-=6,15
 double alpha=2.285; //+-=10
 double H=-0.1270;  //+-=57
 double O=0.4484;  //+-=62

double N=pow(alpha,4+3)/(6*5*4*3*2*1.0);
 double totel=0.0;
int jakovali=128;

char FILENAME[80]="5unitcell.Si";

if(argv<2) 
  {
  printf("No command line parameters given, calculating total density\n");
  n1s=2;n3s=2;n2s=2;n2p=6;n3p=2;
  }
else if(argv==2 && strcmp(argc[1],"3p")==0)
  {
  printf("Calculating 3p density\n");
  n1s=0;n3s=0;n2s=0;n2p=0;n3p=2;
  }
else
  {
  printf("No proper command line parameters given, calculating total density\n");
  printf("Use [3p] parameter to calculate 3p density\n");
  n1s=2;n3s=2;n2s=2;n2p=6;n3p=2;
  }

fp=fopen(FILENAME, "rt");  

//coordinates in Ongstroms

x0=(2*unitcellx);  //upper left corner of the middle unitcell  
y0=(2*unitcelly);
z0=(2*unitcellz);

//delta =unitcelllength/number of dividing units
dx=unitcellx*1/jakovali; 
dy=unitcelly*1/jakovali;
dz=unitcellz*1/jakovali;
npts=0;
nfac[0]=sqrt(ksiexp[0]*ksiexp[0]*ksiexp[0]*4);
nfac3p[0]=sqrt(ksiexp3p[0]*ksiexp3p[0]*ksiexp3p[0]*ksiexp3p[0]*ksiexp3p[0]*4/3);
for(n=1;n<8;n++){
/* Norm factors, Deutsch PrB eq. 3 latter*/
nfac[n]=sqrt(ksiexp[n]*ksiexp[n]*ksiexp[n]*ksiexp[n]*ksiexp[n]*ksiexp[n]*ksiexp[n]*0.17777778);
nfac3p[n]=ksiexp3p[n]*ksiexp3p[n]*ksiexp3p[n]*ksiexp3p[n]*ksiexp3p[n]*ksiexp3p[n]*ksiexp3p[n]*ksiexp3p[n]*ksiexp3p[n];
nfac3p[n]=sqrt(nfac3p[n]*0.012698413);
}

while((ch=fgetc(fp))!=EOF)
{
if(ch==10)
npts++;
}
rewind(fp);
for(n=0;n<npts;n++){
fscanf(fp,"%g %g %g %g", &y[0][n], &y[1][n], &y[2][n],&y[3][n]);
}
fclose(fp);
printf("Data read in (%d points) succesfully\n",npts);
fp=fopen("cdensity.dat","wt");  //file that consists only of the densities (e/Ong^3)
fp2=fopen("cdensreal.dat","wt"); //file that consists of coordinates (Ong) with densities
 
/* loops over x,y and z */
for(nx=0;nx<jakovali;nx++){
xx=x0+nx*dx;
	for(ny=0;ny<jakovali;ny++){
	yy=y0+ny*dy;
		for(nz=0;nz<jakovali;nz++){
		zz=z0+nz*dz;
		ksisum=totsum=0.0;
		for(n=0;n<npts;n++){
		/* loop over atoms */
//Ongstroms -> a0=a_b=Bohr radius (atomic units) 
ax=(y[0][n]-xx)/a_b;
ay=(y[1][n]-yy)/a_b;
az=(y[2][n]-zz)/a_b;

 x2=ax*ax;
 y2=ay*ay;
 z2=az*az;

r=sqrt(x2+y2+z2);
r2=r*r;
r3=r2*r;
r4=r2*r2;

if(r<10){
ksi[0]=nfac[0]*exp(-ksiexp[0]*r*kappa[1])*Y00; /* m=0 1s-orbitaali, -> kappa=1 */
ksi2s[0]=nfac[0]*exp(-ksiexp[0]*r*kappa[0])*Y00;
ksi2p[0]=nfac3p[0]*r*kappa[0]*exp(-ksiexp[0]*r*kappa[0])*Y00;
ksi0[0]=nfac[0]*exp(-ksiexp[0]*r)*Y00;
ksi3p[0]=nfac3p[0]*r*kappa[1]*Y00*exp(-ksiexp3p[0]*r*kappa[1]);
for(m=1;m<8;m++){
/* 1s, 2s, 3s use one ksi-set; 2p, 3p another one */
/* NOTE: 2s, 2p: kappa[0]; 3s, 3p: kappa[1]*/ 
ksi0[m]=nfac[m]*r2*exp(-ksiexp[m]*r)*Y00;
ksi[m]=nfac[m]*r2*kappa[1]*kappa[1]*exp(-ksiexp[m]*r*kappa[1])*Y00;
ksi3p[m]=nfac3p[m]*r3*kappa3[1]*Y00*exp(-ksiexp3p[m]*r*kappa[1]); /* Laske Slaterin orbitaalit */
ksi2s[m]=nfac[m]*r2*kappa[0]*kappa[0]*exp(-ksiexp[m]*r*kappa[0])*Y00;
ksi2p[m]=nfac3p[m]*r3*kappa3[0]*exp(-ksiexp[m]*r*kappa[0])*Y00;
}



phi1=phi2=phi2p=phi3=phi3p=0;
for(m=0;m<8;m++){		/* eq 3, Y_00 inside */
phi3+=ksi[m]*s3[m]; //3s		/* Lin. combination of Slater orbitals */
phi3p+=ksi3p[m]*s3p[m]; //3p
phi1+=ksi0[m]*s1[m];  //1s 
phi2+=ksi2s[m]*s2[m];  //2s
phi2p+=ksi2p[m]*s2p[m]; //2p
}
totsum+=n3s*phi3*phi3*kappa3[1]+n3p*phi3p*phi3p*kappa3[1]; // eq 2, Y_00 inside 
totsum+=n1s*phi1*phi1+n2s*phi2*phi2*kappa3[0]+n2p*phi2p*phi2p*kappa3[0];
ksisum-=(r*exp(-alpha*r)*O*N*2.0*ax*ay*az)*y[3][n]; // -= sign convention,y[3]=bond +/- 
ksisum+=exp(-alpha*r)*V*((x2*x2+y2*y2+z2*z2)-0.6*r4)*H*N;
}
		}
/*ksisum=6.7483426*ksisum;*/ /* 1/a_b^3 */
totsum+=ksisum;
//charge density (e/a0^3) -> (e/Ongs^3)
totsum/=(a_b*a_b*a_b);  
 totel+=totsum*(dx*dy*dz);
fprintf(fp,"%g\n",totsum);
fprintf(fp2,"%g %g %g %g\n",xx,yy,zz,totsum);
		}
		}
	}

fclose(fp);
fclose(fp2);
printf("Total charge=%g\n",totel);

}






