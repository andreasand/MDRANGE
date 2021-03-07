#include<stdio.h>
#include<math.h>

main() {

int n,npts, maxpts=800;
float SP;
char ch;
FILE *fp;
FILE *gp;
char FILENAME[80];
double a_b=0.529177343575629789;  //Ongs
double pi=3.14159265358979323846;
double hartree=27.2116; //eV
printf("\n");
printf("Give the file name: ");
scanf("%s", FILENAME);
fp=fopen(FILENAME,"rt");
gp=fopen("chargedens_rs.in", "wt");
npts=0;
while((ch=fgetc(fp))!=EOF)
{
if(ch==10)
npts++;
}
rewind(fp);
for(n=0;n<npts;n++)
{
fscanf(fp,"%g", &SP);
SP=SP*(a_b*a_b*a_b);	/* e/Å^3 -> e/a_b^3 atomic units */
SP=pow((3.0/(4.0*pi*SP)),0.33333);   /* rho -> r_s, vakio (3/4 pi)^1/3 */
SP=log(SP);  //take logarithm to make the interpolation in MDRANGE more accurate
fprintf(gp,"%g\n",SP);	
}
fclose(fp);
fclose(gp);
}
