/*
        Zbl96 is a program for calculating electronic and nuclear
        stopping powers according to the semiempirical model of
        Ziegler, Biersack and Littmark.
        
        This program is based on the version 96 of Srim-code.
        
        Version 0.9 written by K. Arstila 16.10.1996

        DO NOT DISTRIBUTE OUTSIDE THE ACCELERATOR LABORATORY OF THE
        UNIVERSITY OF HELSINKI WITHOUT PERMISSION OF THE AUTHOR, 
	EXCEPT AS PART OF THE MDRANGE PACKAGE.

                        Kai.Arstila@Helsinki.FI

*/

#include "zbl96.h"

int main(int argc,char *argv[])
{
   double scoef[ROWS][COLS];
   double m1,m2,E,rho,x,min,max,step,S,sunit,xunit;
   int z1,z2;
   unsigned int flag;
   int i;

   readscoef(scoef);

   readparms(argc,argv,&z1,&z2,&m1,&m2,&rho,&min,&max,&step,&flag,scoef);

   sunit = NA*rho/(m2*1.0e25);
   xunit = 1.0;

   switch(flag & SUNIT){
      case EV_A:
         sunit = 100.0*NA*rho/(m2*1.0e25);      
         break;
      case KEV_NM:
         sunit = NA*rho/(m2*1.0e25);
         break;
      case KEV_UM:
         sunit = 1000.0*NA*rho/(m2*1.0e25);
         break;
      case MEV_MM:
         sunit = 1000.0*NA*rho/(m2*1.0e25);
         break;
      case KEV_UG_CM2:
         sunit = NA/(m2*1e24);
         break;
      case MEV_MG_CM2:
         sunit = NA/(m2*1e24);
         break;
      case KEV_MG_CM2:
         sunit = 1000.0*NA/(m2*1e24);
         break;
      case EV_1E15ATOMS_CM2:
         sunit = 1.0;
         break;
   }
   switch(flag & XUNIT){
      case EV:
         xunit = 1000.0;
         break;
      case KEV:
         xunit = 1.0;
         break;
      case MEV:
         xunit = 0.001;
         break;
      case V0:
         xunit = 1.0;
         break;
      case BETA:
         xunit = 0.0072974;
         break;
      case M_S:
         xunit = 2187673.0;
         break;
      case CM_S:
         xunit = 218767300.0;
         break;
   }

   if(flag & DSA){
      for(x=min,i=0;x<=max;x+=step,i++);
      printf("       %i %10.2f\n",i,rho);
   }

   for(x=min;x<=max;x+=step){
      switch(flag & ENERGY){
         case FALSE:
            E = 25.0*x*x/(xunit*xunit);
            break;            
         default:
            E = x/(xunit*m1);
            break;
      }

      if(sqrt(E)/(5.0*137.035) > 1.0)
         fatal_error(8);
      switch(z1){
         case 1:
            S = pstop(z2,E,scoef);
            break;
         case 2:
            S = hestop(z2,E,scoef);
            break;
         default:
            S = histop(z1,z2,E,scoef);
            break;
      }
      switch(flag & NUCLEAR){
         case N_ONLY:
            S = nuclear(z1,z2,m1,m2,E*m1);
            break;
         case N_BOTH:
            S += nuclear(z1,z2,m1,m2,E*m1);
            break;
         case N_NO:
            break;            
      }
      printf("%12.3e %12.3e\n",x,S*sunit);
   }

   exit(0);
}
void readparms(int argc,char *argv[],int *z1,int *z2,double *m1,double *m2,
               double *rho,double *min,double *max,double *step,
               unsigned int *flag,double scoef[][COLS])
{
   int a,i,j=0;
   char s[5];

   *flag = DEFAULT;

   if(argc==1)
      usage();

   for (i=1;i<argc;i++) {
      if (!strcmp(argv[i],"-h") || argc==1)
         usage();

      if (!strcmp(argv[i],"-n")) 
         *flag = N_ONLY | (*flag & ~NUCLEAR);
      if (!strcmp(argv[i],"-nel"))
         *flag = N_BOTH | (*flag & ~NUCLEAR);
      if (!strcmp(argv[i],"-el"))
         *flag = N_NO | (*flag & ~NUCLEAR);

      if (!strcmp(argv[i],"-e"))
         *flag = KEV | (*flag & ~XUNIT);
      if (!strcmp(argv[i],"-v"))
         *flag = V0 | (*flag & ~XUNIT);

      if (!strcmp(argv[i],"-dsa"))
         *flag |= DSA;


      if (!strcmp(argv[i],"-1") || !strcmp(argv[i],"-eV/Å"))
         *flag = EV_A | (*flag & ~SUNIT);
      if (!strcmp(argv[i],"-2") || !strcmp(argv[i],"-keV/nm"))
         *flag = KEV_NM | (*flag & ~SUNIT);
      if (!strcmp(argv[i],"-3") || !strcmp(argv[i],"-keV/µm"))
         *flag = KEV_UM | (*flag & ~SUNIT);
      if (!strcmp(argv[i],"-4") || !strcmp(argv[i],"-MeV/mm"))
         *flag = MEV_MM | (*flag & ~SUNIT);
      if (!strcmp(argv[i],"-5") || !strcmp(argv[i],"-keV/(µg/cm²)"))
         *flag = KEV_UG_CM2 | (*flag & ~SUNIT);
      if (!strcmp(argv[i],"-6") || !strcmp(argv[i],"-MeV/(mg/cm²)"))
         *flag = MEV_MG_CM2 | (*flag & ~SUNIT);
      if (!strcmp(argv[i],"-7") || !strcmp(argv[i],"-keV/(mg/cm²)"))
         *flag = KEV_MG_CM2 | (*flag & ~SUNIT);
      if (!strcmp(argv[i],"-8") || !strcmp(argv[i],"-eV/(1e15 atoms/cm²)"))
         *flag = EV_1E15ATOMS_CM2 | (*flag & ~SUNIT);

      if (!strcmp(argv[i],"-10") || !strcmp(argv[i],"-eV"))
         *flag = EV | (*flag & ~XUNIT);
      if (!strcmp(argv[i],"-11") || !strcmp(argv[i],"-keV"))
         *flag = KEV | (*flag & ~XUNIT);
      if (!strcmp(argv[i],"-12") || !strcmp(argv[i],"-MeV"))
         *flag = MEV | (*flag & ~XUNIT);
      if (!strcmp(argv[i],"-13") || !strcmp(argv[i],"-v0"))
         *flag = V0 | (*flag & ~XUNIT);
      if (!strcmp(argv[i],"-14") || !strcmp(argv[i],"-beta"))
         *flag = BETA | (*flag & ~XUNIT);
      if (!strcmp(argv[i],"-15") || !strcmp(argv[i],"-m/s"))
         *flag = M_S | (*flag & ~XUNIT);
      if (!strcmp(argv[i],"-16") || !strcmp(argv[i],"-cm/s"))
         *flag = CM_S | (*flag & ~XUNIT);

      if(argv[i][0] != '-'){        /* Z1 Z2 min max step -string */
         if(isdigit(argv[i][0])){
            while(isdigit(argv[i][j]) && (argv[i][j])!='\0')
               j++;
            if(argv[i][j] == '\0'){
               *z1 = atof(argv[i]);
               if(*z1 < 1 || *z1 > 92)
                  fatal_error(4);
               *m1 = scoef[*z1][3];
            } else {
               strncpy(s,argv[i],j);
               s[j] = '\0';
               a = atoi(s);
               get_element(argv[i]+j,a,z1,m1,scoef);
            }
         } else 
            get_element(argv[i],MAI,z1,m1,scoef);

         i++;        /* Dirty hack: we increase the counter inside the loop! */

         if(i>=argc)
            fatal_error(1);

         if(isdigit(argv[i][0])) {
            *z2 = atoi(argv[i]);
             if(*z2 < 1 || *z2 > 92)
               fatal_error(5);            
            *m2 = scoef[*z2][4];
         } else
            get_element(argv[i],NATURAL,z2,m2,scoef);            

         i++;

         if(i>=argc)
            fatal_error(1);
         
         *min = atof(argv[i]);
         i++;

         if(i>=argc){
            *max = *min;
            *step = *min;
         } else {         
            *max = atof(argv[i]);
            i++;

            if(i>=argc)
               fatal_error(1);
         
            *step = atof(argv[i]);
         }
      }
   }
   if(*z1 < 1 || *z1 > 92)
      fatal_error(4);
   if(*z2 < 1 || *z2 > 92)
      fatal_error(5);
   if(*min < 0 || *max < 0)
      fatal_error(3);
   if(*step <= 0.0)
      fatal_error(7);          
   if(*max < *min)
      fatal_error(2);
   
   *rho = scoef[*z2][5];
}
void fatal_error(int errno)
{
   fprintf(stderr,"      Error: %s\n\n",err_strings[errno]);
   
   exit(errno);
}
void usage(void)
{
   fprintf(stderr,"zbl96 ver 0.9 © K. Arstila and J.F. Ziegler\n");
   fprintf(stderr,"Usage: zbl96 [options] Z1 Z2 min [max step]\n\n");
   fprintf(stderr,"     Z1 and Z2 can be given as chemical symbols\n");
   fprintf(stderr,"     Z1 can be preceded by a number of nucleons\n\n");
   fprintf(stderr,"Options: \n");
   fprintf(stderr,"    -h                        help\n");
   fprintf(stderr,"    -n                        only nuclear stopping\n");
   fprintf(stderr,"    -nel                      both nuclear and electronic stopping\n");
   fprintf(stderr,"    -el                       only electronic stopping (default)\n");
   fprintf(stderr,"    -e                        stopping as a function of energy\n");
   fprintf(stderr,"    -v                        stopping as a function of velocity (default)\n");
   fprintf(stderr,"    -dsa                      includes number of points and density to output\n");
   fprintf(stderr,"    -1 -eV/Å                  stopping unit is eV/Å\n");
   fprintf(stderr,"    -2 -keV/nm                stopping unit is keV/nm (default)\n");
   fprintf(stderr,"    -3 -keV/µm                stopping unit is keV/µm\n");
   fprintf(stderr,"    -4 -MeV/mm                stopping unit is MeV/mm\n");
   fprintf(stderr,"    -5 -'keV/(µg/cm²)'        stopping unit is keV/(µg/cm²)\n");
   fprintf(stderr,"    -6 -'MeV/(mg/cm²)'        stopping unit is MeV/(mg/cm²)\n");
   fprintf(stderr,"    -7 -'keV/(mg/cm²)'        stopping unit is keV/(mg/cm²)\n");
   fprintf(stderr,"    -8 -'eV/(1e15 atoms/cm²)' stopping unit eV/(1e15 atoms/cm²)\n");
   fprintf(stderr,"    -10 -eV                   energy unit is eV\n");
   fprintf(stderr,"    -11 -keV                  energy unit is keV (default)\n");
   fprintf(stderr,"    -12 -MeV                  energy unit is MeV\n");
   fprintf(stderr,"    -13 -v0                   velocity unit is Bohr velocity (default)\n");
   fprintf(stderr,"    -14 -beta                 velocity unit is relative to the velocity of light\n");
   fprintf(stderr,"    -15 -m/s                  velocity unit is m/s\n");
   fprintf(stderr,"    -16 -cm/s                 velocity unit is cm/s\n");

   exit(0);
}
void get_element(char *s,int a,int *z,double *m,double scoef[][COLS])
{
   FILE *fp;
   char S[3];
   int N,Z,A,len,found = FALSE;
   double M;
   
   len = strlen(s);
   
   fp = fopen(F_ELEMENTS,"r");
  
   while(!found && fscanf(fp,"%i %i %i %s %lf",&N,&Z,&A,S,&M)==5){
      if(strncmp(S,s,len) == 0 && (a == A || a == MAI || a == NATURAL))
         found = TRUE;
   }
   fclose(fp);

   if(!found){
      switch(a){
         case MAI:
            fatal_error(5);
            break;
         case NATURAL:
            fatal_error(5);         
            break;
         default:
            fatal_error(6);
            break;
      }
   }
    
   *z = Z;
   switch(a){
      case MAI:
         if(*z < 1 || *z > 92)
            fatal_error(5);
         *m = scoef[*z][3];
         break;
      case NATURAL:
         if(*z < 1 || *z > 92)
            fatal_error(5);
         *m = scoef[*z][4];
         break;
      default:
         *m = M/1000000.0;
         break;
   }
}
