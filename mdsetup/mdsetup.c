/*
  mdsetup - Setup program for the input to the range calculation code mde

  © 1996- J. K. Sillanpää and Kai Nordlund (kai.nordlund@helsinki.fi)

  Note: mdsetup is sort of beta still, it has not been thoroughly
  tested. If you find bugs in it, please let us know. In the meantime,
  you can always circumvent mdsetup, and edit the mdh input files 
  by hand.

*/

#define VERSION   "V0.96"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"MDpt.h"
#include"local.h"
#include<ctype.h>
#include<math.h>
#include<time.h>
#include<unistd.h>

#define True 1
#define False 0

typedef int logical;

main() {
  int lkm, n, j, k, l, p, h, s=0,t, m, mo, grid, lkma, mmax, comp, multi=0,
    seed, layers, rtype=1, atrpot=0, donetypes=0, szbl=0, expert=0, Ncalc=10000, 
    Ninistate=10, Trange=0, pakkot=0, help=48, nice=8, maxtypes=0, ilk;
  int  types[1000], reptyp[8];
  int cells[3]={2, 2, 2};
  int Z[9];
  float atoms [2000][4], laydim[2][4]; 

  float OTHER [12] [3], xred [1000] [3];
  const float BCC [2] [3] = {{0, 0, 0}, {0.5, 0.5, 0.5}};
  const float FCC [4] [3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const float DIA [8] [3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5},
			     {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  const float HCP [4] [3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.0, 0.33333333, 0.5}, 
			     {0.5, 0.833333, 0.5}}; 

  float elstop[50] [2];
  float lim[3]={0.5, 0.5, 0.5};
  float reppot[6] [2]= {{0.1818, -3.2}, {0.5099, -0.9423},
			{0.2802, -0.4029}, {0.02817, -0.2016}, {0, 0}, {0, 0}};

  int laycomp[10]; 
  float min,max,dh,rat;
  double boxx,boxy,boxz;
  int nums[10][5],numsa[10],lukum,nim,suu,pie, ntyp [10] ;
  char fname[140];
  float sput_surf=0.0,Dose=1e10;
  int Damage=0,sputtering=0,scaling=0,bscaling=0;
  int inpostyp;
  float multi_a[10],multi_b[10],multi_c[10]; /*multilayer unitcell dims..*/
  int multi_Z[10][5];  /*10 layers,5-compound max*/
  float speed, force, transl[3], totalth, polysize=0, 
    polysizesigma, delta, E0, Theta0=0, minz=-3, maxz=-3,
    Thetamax=0,Fii0=0, Fiimax=0, miny=0.5, minx=0.5, 
    maxx=1.5,maxy=1.5, Tini=300, 
    Timeini=300, elscale=1, repscale=1, R0rec=3, Rmrec=4, Rmini=5.33, R0ini=4.3; 

  double pi=3.1415926535;
  float adens;
  double a,b,c, incrspeed, topspeed;
  double mass[9], dtemp[4]={0, 0, 0, 0};
  FILE *fp;
  FILE *gp;
  FILE *ep;
  FILE *hp;
  FILE *cp;
  FILE *xp;

  char unit2[3]="-8";
  char file[80]="coords.in";
  char file1[80]="coords.in";
  char file2[80]="elstop.in";
  char file3[80]="param.in";
  char file4[80]="reppot.0.1.in";
  char file5[80]="runMDH";
  char OK[2];
  char ch, projectile[4], lattice[5],zsym[3], buf[80], speeddat[80];
  char target[5], nicei[4], name[5], numb[3];

  /* if you have a different version of ZBL, */
  /*  modify these lines accordingly */

  char inizbl[100] = BINDIR"zbl96 ";    
  char zbl[200];			      
  char nices[90]="nice -";
  char chmod[90]="chmod 777 ";
  char mdh[200]="mdh ";
  char end[] = " >temp.in";
  char compoundname [120];
  char compounds [120];
  char endcomp[]="compounds.txt";
  time_t walltime;

  logical elstoperror=False;

  void givehelp(int help);

  time(&walltime);

  strcat(compounds,DATADIR);
  strcat(compounds,endcomp);
  printf("\n");
  printf("MDsetup %s\n\n",VERSION);
  printf("Are you an EXPERT USER ?\n");
  printf("Declaring yourself as one will result in a multitude of technical\n");
  printf("questions, which will otherwise be bypassed by default values !\n"); 
  printf("(Y/N) [N]: ");
  if(fgets(buf,79,stdin)!=NULL)
    sscanf(buf,"%s",OK);
  if(OK[0]=='Y' || OK[0]=='y')
    expert=1;
	
  printf("\n");
  printf("#----------------------------------------------------------------#\n");
  printf("# 			The Target Material                      #\n"); 
  printf("#----------------------------------------------------------------#\n");
  printf("\n");
  /*
    Because the simulation of multilayered targets is somewhat buggy,
    these features have been commented out. If you want to simulate such 
    targets, write to us for modified simulation code.
    jksillan@vesuri.helsinki.fi
    kai.nordlund@helsinki.fi*/
  OK[0]='A';
  printf("Is the target material multilayered (Y/N) [N]: ");
  if(fgets(buf,79,stdin)!=NULL)
    sscanf(buf,"%c",OK);
  if(OK[0]=='Y' || OK[0]=='y') {
    laydim[0][0]=-1e5;
    multi=1;
    printf("How many layers are there (<=4): ");
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf,"%d",&layers);    
  }
  else
   
    layers=1;
     

  for(t=0;t<layers;t++) {
    comp=1;
    printf("\nIs the %d. target material an element or a compound? \n", t+1); 
    printf("1 element 2 compound in the list 3 other compound [1]: ");
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf,"%d",&comp);
    laycomp[t]=comp; printf("%d. layer, comp=%d\n",t+1,laycomp[t]);



    if(comp==1) {
      for(;;) {
	l=0;
	printf("\nPlease enter the atomic number or chemical symbol of the lattice atoms: "); 
	scanf("%s",zsym);
	zsym[0]=toupper(zsym[0]);
	zsym[1]=tolower(zsym[1]);
	if(48 <= zsym[0] && zsym[0] <=57) {
	  if(zsym[1]==10 || zsym[1]==0) {
	    Z[1+donetypes]=zsym[0]-48;
	    maxtypes=lkma=l=1;
	  }
	  else {
	    Z[1+donetypes]=zsym[1]-48+((zsym[0]-48)*10);
	    maxtypes=lkma=l=1;
	  }
	}
	else {
	  if(zsym[1]==10 || zsym[1]==0) {
	    for(n=0;n<108;n++) {
	      strcpy(name,p_table[n].symbol);
	      if(name[0]==zsym[0] && name[1] < 97) {
		maxtypes=l=lkma=1;
		Z[1+donetypes]=n+1;
	      }
	      if(l==1)
		break;
	    }
	  }
	  else {
	    for(n=0;n<108;n++) {
	      strcpy(name,p_table[n].symbol);
	      if(name[0]==zsym[0] && name[1]==zsym[1]) {
		maxtypes=l=lkma=1;
		Z[1+donetypes]=n+1;
	      }
	      if(l==1)
		break;
	    }
	  }
	  if(l==0)
	    printf("Element %s not found!\n", zsym);
	}
	if(l==1) {
	  mass[1+donetypes]=p_table[Z[1+donetypes]-1].atomic_mass;
	  if(mass[1+donetypes]<0)
	    mass[1+donetypes]=-mass[1+donetypes];
	  strcpy(name,p_table[Z[1+donetypes]-1].symbol);
	  printf("%s , Z=%d\n", name, Z[1+donetypes]);

	  multi_Z[t][0]=Z[1+donetypes]; 
	       

	  OK[0]='A';
	  types[0]=1+donetypes; 
	  if(fgets(buf,79,stdin)!=NULL) {
	  }
	  printf("Is this OK (Y/N) [Y] :");
	  if(fgets(buf,79,stdin)!=NULL)
	    sscanf(buf,"%c",OK);
	  if(OK[0]=='N' || OK[0]=='n')
	    l=0;
	}
	if(l==1)
	  break;
      }

      strcpy(lattice,p_table[Z[1+donetypes]-1].crystal);
      a=p_table[Z[1+donetypes]-1].a;
      l=0;
      sprintf(target, " %d ", Z[1+donetypes]);

      for (;;) {
	OK[0]='A';
	printf("\nAtomic mass: %g\n",mass[1+donetypes]);
	printf("Crystal type: %s\n",lattice);
	printf("Length of unit cell in x dimension: %g\n",a);
	OK[0]='A';
	if(lattice[0]!='B'&&lattice[0]!='F'&&lattice[0]!='D'&&lattice[0]!='H'&&lattice[0]!='S') {
	  printf("Warning! Mdsetup DOES NOT CURRENTLY SUPPORT this crystal type.\n");
	  printf("Unless the crystal type is changed, \n");
	  printf("You will have to make the coordinates file yourself\n");
	}
	printf("Is this OK (Y/N) [Y]: ");
	if(fgets(buf,79,stdin)!=NULL);
	sscanf(buf, "%c",OK);
	if(lattice[0]=='B') grid=1;
	if(lattice[0]=='F') grid=2;
	if(lattice[0]=='D') grid=3;
	if(lattice[0]=='H') grid=4;
	if(lattice[0]=='S') grid=5;

	if(OK[0]!='N' && OK[0]!='n') {
	  l=1; 
	}
	else {
	  printf("\n");
	  printf("Please enter the type of the lattice: \n");
	  printf("1) BCC \n");
	  printf("2) FCC \n");
	  printf("3) DIA \n");
	  printf("4) HCP \n");
	  printf("5) SC \n");
	  printf(": ");
	  scanf("%d", &grid);
	  if(grid==1)
	    strcpy(lattice,p_table[2].crystal);
	  if(grid==2)
	    strcpy(lattice,p_table[9].crystal);
	  if(grid==3)
	    strcpy(lattice,p_table[5].crystal);
	  if(grid==4)
	    strcpy(lattice,p_table[0].crystal);
	  if(grid==5)
	    strcpy(lattice,p_table[83].crystal);
	  printf("Please enter the mass of the lattice atoms(amu): ");
	  scanf("%lf", &mass[1+donetypes]);
	  printf("Please enter the length of the unit cell(Å): ");
	  scanf("%lg", &a);
	  if(fgets(buf,79,stdin)!=NULL) {
	  }
	}
	if (l==1)
	  break;
      }
      b=a;
      if(grid==4) {
	b=sqrt(3.0)*a;
	c=1.6*a;               
	printf("Creating HCP with 4 atom rectangular unit cell\n");
	printf("Please enter height of unit cell(Å) [%lg]: ",c);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%lg", &c);
	printf("OK, rectangular HCP unit cell size is %g %g %g\n",a,b,c);
      }
      else c=a;	 

      dtemp[t]=p_table[Z[1+donetypes]-1].tdebye;
      if (dtemp[t]<0) dtemp[t]=-dtemp[t];
	 
    }
    /* Handle compound selection */

    if(comp==2) {
      cp=fopen(compounds,"rt");
      if (cp==NULL) {
	printf("Failed to open compounds file %s\n",compounds);
      }
      printf("\n");
      m=mmax=1;
      for(;;) {
	while((ch=fgetc(cp))!= '%' && (m%40)!=0) {
	  while((ch=fgetc(cp)) != '$') { }
	  fscanf(cp,"%s", compoundname);
	  printf("%d %s", m, compoundname); fflush(stdout);
	  if((m%2)==0)
	    printf("\n"); 
	  else
	    printf("\t");
	  m++;
	  mmax++;
	  while((ch=fgetc(cp)) != 35) {
	  }
	}
	printf("\n");
	mmax--;
	printf("Which material [more compounds]: ");
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%d", &m);
	printf("\n");
	if(m>0 && m<=mmax)
	  break;
	if((mmax%40)!=0) {
	  m=mmax=1;
	  rewind(cp);
	}
 
      }
      rewind(cp);
 
      for(n=0;n<m;n++) {
	while((ch=fgetc(cp)) != 36) {
	}   
      }
      fscanf(cp,"%s", compoundname);
      fscanf(cp,"%s", lattice);
      if(lattice[0]=='B') { /* BCC */
	grid=1;
	lkma=2;
      } 
      if(lattice[0]=='F') { /* FCC */
	grid=2;
	lkma=4;
      }
      if(lattice[0]=='D') { /* Diamond */
	grid=3;
	lkma=8;
      }
      if(lattice[0]=='H') { /* HCP */
	grid=4;
	lkma=4;
      }
      if(lattice[0]=='O') {  /* OTHER unknown type */
	grid=6;
      }
      fscanf(cp,"%lf", &dtemp[t]);
      fscanf(cp,"%d",&maxtypes);
      fscanf(cp,"%lf",&a);
      printf("According to the database, material has lattice type %s\n",lattice);
      printf("Please enter the lattice length [%lf]: ", a);
      if(fgets(buf,79,stdin)!=NULL) sscanf(buf,"%lf",&a);

      if(grid==4 || grid==6) {
	fscanf(cp,"%lf",&c);
	printf("Please enter the lattice height [%lf]: ", c);
	if (fgets(buf,79,stdin)!=NULL) sscanf(buf,"%lf",&c);

	if (grid==4) {
	  printf("Unit cell size for type HCP is %lf %lf %lf and N %d\n",a,b,c,lkma);
	}
	if(grid==6) {
	  fscanf(cp,"%lf",&b);
	  printf("Please enter the lattice depth [%lf]: ", b);
	  if(fgets(buf,79,stdin)!=NULL) sscanf(buf,"%lf",&b);
	  fscanf(cp,"%d",&lkma);
	  printf("Unit cell size for type OTHER is %lf %lf %lf and N %d\n",a,b,c,lkma);
	}
      }
      else { 
	c=a;
	b=a;
      }


      for(n=0;n<lkma;n++) {
	fscanf(cp,"%d",&(types[n]));
	types[n]=types[n]+donetypes;
      }

      for(m=donetypes+1;m<=maxtypes+donetypes;m++)
	{fscanf(cp,"%d",&(Z[m]));multi_Z[t][m-donetypes-1]=Z[m];}

      nums[t][0]=0; for(n=0;n<lkma;n++) if(Z[types[n]]==multi_Z[t][0]) nums[t][0]++;
      nums[t][1]=0; for(n=0;n<lkma;n++) if(Z[types[n]]==multi_Z[t][1]) nums[t][1]++;
      numsa[t]=lkma;


      for(n=1;n<=maxtypes;n++)
	fscanf(cp,"%lf",&(mass[n+donetypes]));
      if(grid==6) {
	for(n=0;n<lkma;n++) {
	  fscanf(cp,"%g", &(OTHER[n][0]));
	  fscanf(cp,"%g", &(OTHER[n][1]));
	  fscanf(cp,"%g", &(OTHER[n][2]));	
	}
      }
      fclose(cp);
      m=mmax=0;

    }



    if(comp==3) {
      /*local dev by JPC never tested on multilayers */
      printf("How many atom types (max=6) \n ");
      scanf("%d",&maxtypes);    
      printf("maxtypes %d \n",maxtypes);
	      
      for(ilk=1;ilk<=maxtypes;ilk++) {   /* loop on atom types */
	for(;;) {
	  l=0;
	  printf("\nPlease enter the atomic number or chemical symbol of the atom type %d: ",ilk); 
	  scanf("%s",zsym);
	  zsym[0]=toupper(zsym[0]);
	  zsym[1]=tolower(zsym[1]);
	  if(48 <= zsym[0] && zsym[0] <=57) {
	    if(zsym[1]==10 || zsym[1]==0) {
	      Z[ilk+donetypes]=zsym[0]-48;
	      l=1;
	    }
	    else {
	      Z[ilk+donetypes]=zsym[1]-48+((zsym[0]-48)*10);
	      l=1;
	    }
	  }
	  else {
	    if(zsym[1]==10 || zsym[1]==0) {
	      for(n=0;n<108;n++) {
		strcpy(name,p_table[n].symbol);
		if(name[0]==zsym[0] && name[1] < 97) {
		  l=1;
		  Z[ilk+donetypes]=n+1;
		}
		if(l==1)
		  break;
	      }
	    }
	    else {
	      for(n=0;n<108;n++) {
		strcpy(name,p_table[n].symbol);
		if(name[0]==zsym[0] && name[1]==zsym[1]) {
		  l=1;
		  Z[ilk+donetypes]=n+1;
		}
		if(l==1)
		  break;
	      }
	    }
	    if(l==0)
	      printf("Element %s not found!\n", zsym);
	  }
	  if(l==1) {
	    mass[ilk+donetypes]=p_table[Z[ilk+donetypes]-1].atomic_mass;
	    if(mass[ilk+donetypes]<0)
	      mass[ilk+donetypes]=-mass[ilk+donetypes];
	    strcpy(name,p_table[Z[ilk+donetypes]-1].symbol);
	    printf("%s , Z=%d\n", name, Z[ilk+donetypes]);
	    printf("Atomic mass: %g\n",mass[ilk+donetypes]);
	    multi_Z[t][ilk-1]=Z[ilk+donetypes]; 
	       

	  }
	  if(l==1)
	    break;
	}
      }


      printf("How many atoms (max 100) ? \n ");
      scanf("%d",&lkma);    

      inpostyp=1 ;

      printf("positions and types entered with keyboard (1) or in file atoms.in (2) :");
      scanf("%d",&inpostyp);
      if(inpostyp==1) {
	for(n=0;n<lkma;n++) {
	  printf("xred yred zred type of atom %d \n ",n+1);
	  
	  scanf("%f%f%f%d",&xred[n][0],&xred[n][1],&xred[n][2],&types[n]);   
	}
	
      }
      else {
	xp = fopen("atoms.in", "rt");
	for(n=0;n<lkma;n++) {
	  fscanf(xp,"%f%f%f%d",&xred[n][0],&xred[n][1],&xred[n][2],&types[n]);   
	}
      }
      
    	

      for(ilk=1;ilk<=maxtypes;ilk++) {   /* loop on atom types */
	nums[t][ilk-1]=0; for(n=0;n<lkma;n++) if(Z[types[n]]==multi_Z[t][ilk-1]) nums[t][ilk-1]++;
	numsa[t]=lkma;
      }
      
      ntyp[t]= maxtypes;

      printf("\n cell supposed tetragonal \n");

      printf(" Please enter the lattice lengths a b c: \n");
      scanf("%lf%lf%lf",&a,&b,&c);    

      /*	printf("  a= %lf \n",a);
	printf(" Please enter the lattice length of b: \n");
	scanf("lf",&b);    
	printf("  b= %lf \n",b);
	printf(" Please enter the lattice length of c: \n");
	scanf("%lf",&c);    
	printf("  c= %lf \n",c);
      */
      printf("Unit cell size  is %lf %lf %lf and N %d\n",a,b,c,lkma);

      grid=7;
    } /* fin comp=3*/




 


    if(multi) {
      printf("\nPlease enter the thickness of the layer (Å) ");
      printf("\n(Note! The program uses multiples of unit cell z-dimension) :");
      scanf("%f", &dh);
      if(t>=0) {
	laydim[1][t]=totalth+dh;
	laydim[0][t]=totalth;
      }
      totalth+=dh;
      if(fgets(buf,79,stdin)!=NULL) {
      }
    }

    if(expert==1 && t==0) {   /*ask only during the first layer query*/
      for(;;) {
	printf("\n");
	printf("    (NEIGBOUR LIST CUT-OFF > POTENTIAL CUT-OFF)\n");
	printf("    Please enter the NEIGHBOUR LIST construction cut-off radius\n");
	printf("    to be used in RECOIL CALCULATION(Å) [%g]: ", Rmrec);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &Rmrec);
	printf("    Please enter the POTENTIAL cut-off radius\n");
	printf("    to be used in RECOIL CALCULATION (Å) [%g]: ", R0rec);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &R0rec);
	if(Rmrec >= R0rec)
	  break;
      }
    }
    donetypes+=maxtypes;

    if(1<=grid && grid<=7) {
      OK[0]='Y';
      printf("\n");
      printf("#--------------------------------------------------------------#\n");
      printf("#		 The Coordinates File			        #\n");
      printf("#--------------------------------------------------------------#\n");
      printf("\n");
      if(expert==1) {
	printf("\nDo you want to create a coordinates file (Y/N) [Y]: ");
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%c",OK);
      }

      lkm=250;


      /* 
	 Note: for crystal symmetry reasons lim[] always has to be an integer,
      */
      if(a<=2*R0rec) { lim[0]=1.0; cells[0]=3; minx=1.0; maxx=2.0; }
      if(a<=R0rec) { lim[0]=2.0; cells[0]=6; minx=2.5; maxx=3.5; }

      if(b<=2*R0rec) { lim[1]=1; cells[1]=3; miny=1.0; maxy=2.0; } 
      if(b<=R0rec) { lim[1]=2.0; cells[1]=6; miny=2.5; maxy=3.5; }

      if(c<=2*R0rec) { lim[2]=1; cells[2]=3; }
      if(c<=R0rec) { lim[2]=2.0; cells[2]=6; }


      multi_a[t]=a;multi_b[t]=b;multi_c[t]=c;  /*unitcell dims for multilayers*/
	    

      if(OK[0]!='N' && OK[0] !='n') {
	printf("\n");
	if(expert==1 && t==0) {  
	  printf("The size of the simulation cell: \n");
	  if(grid==6) {
	    printf("\n");
	    printf("The compound you have chosen has a non-standard lattice type.\n");
	    printf("The unit cell used for this compound consists of %d atoms.\n\n", lkma);
	  }
	       
	  printf("Please enter the number of unit cells in the X-direction [%d]: ", cells[0]); 
	  if(fgets(buf,79,stdin)!=NULL)
	    sscanf(buf,"%d", &cells[0]);

	  printf("Please enter the number of unit cells in the Y-direction [%d]: ", cells[1]); 
	  if(fgets(buf,79,stdin)!=NULL) 
	    sscanf(buf,"%d", &cells[1]);
	  printf("Please enter the number of unit cells in the Z-direction [%d]: ",cells[2]); 
	  if(fgets(buf,79,stdin)!=NULL) 
	    sscanf(buf,"%d", &cells[2]);
	       
	}
	if(t==0) {boxx=a*cells[0];boxy=b*cells[1];boxz=c*cells[2];} /*size of the sim. box*/
	    
	/*if first layer size (=sim box size) is greater than the size for other layers, 
	  this code increases the size of those layers.*/
	if(a*cells[0]<(boxx-0.001)) {printf("Making this layer larger\n");cells[0]=ceil(boxx/a);}
	if(b*cells[1]<(boxy-0.001)) {printf("Making this layer larger\n");cells[1]=ceil(boxy/b);}
	if(c*cells[2]<(boxz-0.001)) {printf("Making this layer larger\n");cells[2]=ceil(boxz/c);}
	    
	printf("\n");
	printf("Creating the lattice, %d * %d * %d unit cells \n", cells[0], cells[1],cells[2]); 
	transl[0]=transl[1]=transl[2]=a/8;
	nim=0;
	if(grid==1) {
	  lkm = 2 * cells[0]*cells[1]*cells[2];
	  transl[0]=transl[1]=transl[2]=a/4;
	  for (n=0; n<cells[0]; n++) {
	    for(k=0; k<cells[1]; k++) {
	      for(j=0; j<cells[2]; j++) {
		for(m=0; m<2; m++) {
		  atoms [s] [0] = transl[0] + (a * (BCC [m] [0])) + a * n;
		  atoms [s] [1] = transl[1] + (a * (BCC [m] [1])) + a * k;
		  atoms [s] [2] = transl[2] + (a * (BCC [m] [2])) + a * j;
		  mo=(s%lkma)+nim;
		  atoms [s] [3] = types[mo];
		  s++;
		}			
	      }
	    }

	  }
	}
	if(grid==3) {
	  lkm=8*cells[0]*cells[1]*cells[2];
	  transl[0]=transl[1]=transl[2]=a/8;
	  for (n=0; n<cells[0]; n++) {
	    for(k=0; k<cells[1]; k++) {
	      for(j=0; j<cells[2]; j++) {
		for(m=0; m<8; m++) {
		  atoms [s] [0] = transl[0]+ (a * (DIA [m] [0])) + a * n;
		  atoms [s] [1] = transl[1]+ (a * (DIA [m] [1])) + a * k;
		  atoms [s] [2] = transl[2]+ (a * (DIA [m] [2])) + a * j;
		  mo=(s%lkma)+nim;
		  atoms [s] [3] = types[mo];
		  s++;
		}			
	      }
	    }

	  }
	}
	if(grid==2) {
	  lkm = 4 * cells[0] * cells[1] * cells[2];
	  transl[0]=transl[1]=transl[2]=a/4;
	  for (n=0; n<cells[0]; n++)  {
	    for(k=0; k<cells[1]; k++) {
	      for(j=0; j<cells[2]; j++) {
		for(m=0; m<4; m++) {
		  atoms [s] [0] = transl[0]+ (a * (FCC [m] [0])) + a * n;
		  atoms [s] [1] = transl[1]+ (a * (FCC [m] [1])) + a * k;
		  atoms [s] [2] = transl[2]+ (a * (FCC [m] [2])) + a * j;
		  mo=(s%lkma)+nim;
		  atoms [s] [3] = types[mo];
		  s++;
		}
	      }
	    }

	  }
	}

	if(grid==4) {
	  lkm = 4 * cells[0] * cells[1] * cells[2];
	  transl[0]=a/4;
	  transl[1]=b/12;
	  transl[2]=c/4;

	  for (n=0; n<cells[0]; n++) {
	    for(k=0; k<cells[1]; k++) {
	      for(j=0; j<cells[2]; j++) {
		for(m=0; m<4; m++) {
		  atoms [s] [0] = transl[0] + (a * (HCP [m] [0])) + a * n;
		  atoms [s] [1] = transl[1] + (b * (HCP [m] [1])) + b * k;
		  atoms [s] [2] = transl[2] + (c * (HCP [m] [2])) + c * j;
		  mo=(s%lkma)+nim;
		  atoms [s] [3] = types[mo];
		  s++;
		}
	      }
	    }

	  }
	}

	if(grid==5) {
	  lkm=cells[0]*cells[1]*cells[2];
	  transl[0]=transl[1]=transl[2]=a/2;
	  for(n=0;n<cells[0];n++) {
	    for(k=0;k<cells[1];k++) {
	      for(j=0; j<cells[2]; j++) {
		atoms [s] [0]= a*n  + transl[0];
		atoms [s] [1]= b*k + transl[1];
		atoms [s] [2]= c*j  + transl[2];
		atoms [s] [3] = 1+nim; 
		s++;
	      }
	    }
	  }
	}


	if(grid==6) {
	  lkm=lkma*cells[0]*cells[2]*cells[1];
	  for (n=0; n<cells[0]; n++) {
	    for(k=0; k<cells[1]; k++) {
	      for(j=0; j<cells[2]; j++) {
		for(m=0; m<lkma; m++) {
		  atoms[s][0]=transl[0]+(a*OTHER[m][0]) + a*n;
		  atoms[s][1]=transl[1]+(b*OTHER[m][1])+  b*k;
		  atoms[s][2]=transl[2]+(c*OTHER[m][2]) + c*j;
		  atoms[s][3]=types[m+nim];
		  s++; 
		}
	      }
	    }
	  }
	}

	if(grid==7) {
	  transl[0]=a/20 ;
	  transl[1]=b/20 ;
	  transl[2]=c/20 ;
	  lkm=lkma*cells[0]*cells[2]*cells[1];
	  for (n=0; n<cells[0]; n++) {
	    for(k=0; k<cells[1]; k++) {
	      for(j=0; j<cells[2]; j++) {
		for(m=0; m<lkma; m++) {
		  atoms[s][0]=transl[0]+(a*xred[m][0]) + a*n;
		  atoms[s][1]=transl[1]+(b*xred[m][1])+  b*k;
		  atoms[s][2]=transl[2]+(c*xred[m][2]) + c*j;
		  atoms[s][3]=types[m];
		  s++; 
		}
	      }
	    }
	  }
	}



	  printf("\n");
	  if(multi==1) {
	    sprintf(numb,".%d",t+1);
	    strcat(file1,numb);
	  }
	  if((fp=fopen(file1, "rt"))!=NULL) {
	    printf("The default coordinates file (%s) already exists!\n", file1); 
	    printf("Where do you want to save the coordinates [%s]: ",file1); 
	    if(fgets(buf,79,stdin)!=NULL)
	      sscanf(buf,"%s", file1);
	    printf("\n");
	    fclose(fp);
	  }

	  if(multi) 
	    {
	      printf("Scaling set for layer boundaries\n");
	      scaling=1;
	    }
	  if(expert==1 && t==0)
	    {
	      if(multi) 
		{
		  scaling=1;
		}
	      scaling=1;
	    }
	  fp = fopen(file1, "wt");

	  {
	    fprintf(fp,"V %g %g %g %d %d %d\n",a,b,c,cells[0],cells[1],cells[2]); 
	    /*REAL size for this cell, does not have to be in, if basic mdrange stuf is used*/
	    fprintf(fp,"D %g \n",0.0);  /*dose in this cell used only when Damage modelling is on*/
	    /*if the bscaling is used to change the cell size, coordinates have to be origo-centered!!!*/
	  }
	    
	  for(h=0;h<lkm;h++) {
	    fprintf(fp, "%g %g %g %g \n", atoms [h] [0], atoms [h] [1], atoms [h] [2], atoms [h] [3]); 
	  }
	  fclose(fp);
	
	
      }
      strcpy(file1,file);
      s=0;
    }
  }
  printf("\n");
  printf("#----------------------------------------------------------------#\n");
  printf("#               The Implanted Ion & Its Interactions             #\n");
  printf("#----------------------------------------------------------------#\n");
  for(;;) {
    OK[0]='A';
    l=0;
    printf("\nPlease enter the atomic number or chemical symbol of the projectile: "); 
    scanf("%s", zsym);
    zsym[0]=toupper(zsym[0]);
    zsym[1]=tolower(zsym[1]);
    if(48 <= zsym[0] && zsym[0] <=57) {
      if(zsym[1]==10 || zsym[1]==0) {
	Z[0]=zsym[0]-48;
	l=1;
      }
      else {
	Z[0]=zsym[1]-48+((Z[0]=zsym[0]-48)*10);
	l=1;
      }
    }
    else {
      if(zsym[1]==10 || zsym[1]==0) {
	for(n=0;n<108;n++) {
	  strcpy(name,p_table[n].symbol);
	  if(name[0]==zsym[0] && name[1]<97) {
	    Z[0]=n+1;
	    l=1;
	  }
	  if(l==1)
	    break;  
	}
      }
      else {
	for(n=0;n<108;n++) {
	  strcpy(name,p_table[n].symbol);
	  if(name[0]==zsym[0] && name[1]==zsym[1]) {
	    l=1;
	    Z[0]=n+1;
	  }
	  if(l==1)   
	    break;
	}
      }
      if(l==0)
	printf("Element %s not found!\n", zsym);
    }
    if(l==1)
      break;
  }
        
  mass[0]=p_table[Z[0]-1].abundant_mass;
  if(mass[0]==0) {
    printf("Please enter the mass of the projectile (amu): ");
    scanf("%lf", &mass[0]);
  }
  else  {
    printf("Retrieving the mass of the most stable isotope of %s \n",p_table[Z[0]-1].symbol);
  }
  if(fgets(buf,79,stdin)!=NULL) { }

  printf("%s, Z=%d, mass=%g\n", p_table[Z[0]-1].symbol, Z[0], mass[0]);
  printf("Is this OK (Y/N) [Y] :");
  if(fgets(buf,79,stdin)!=NULL)
    sscanf(buf,"%c",OK);
  if(OK[0]=='N' || OK[0]=='n') {
    for(;;)/*projectile check-up*/ {
      printf("\n");
      printf("Please enter the atomic number of the projectile: ");
      scanf("%d",&Z[0]);
      printf("Please enter the mass of the projectile (amu): ");
      scanf("%lf", &mass[0]);
      printf("\n");
      if(fgets(buf,79,stdin)!=NULL) {
      }
      OK[0]='A';
      printf("%s, Z=%d, mass=%g\n", p_table[Z[0]-1].symbol, Z[0], mass[0]);
      printf("Is this OK (Y/N) [Y] :");
      if(fgets(buf,79,stdin)!=NULL)
	sscanf(buf,"%c",OK);
      if(OK[0]!='N' && OK[0]!='n')
	break;
    }
  }
	



  sprintf(projectile, "%d", Z[0]);
  printf("\nPlease enter the recoil energy (keV): ");
  fgets(buf,90,stdin); sscanf(buf,"%g", &E0);

  if (E0<=0.0) { 
    printf("I can not handle E <= 0. Using E = 0.001 keV\n");
    E0=0.001;
  }

  topspeed=(E0*1.07e-6/mass[0])+1;
  topspeed = 1.1*sqrt(1-pow(topspeed,(double) -2.0))*150;
  incrspeed=topspeed/50;
  E0 = E0 * 1000;
  sprintf(speeddat, " 0.0 %lg %lg \0", topspeed, incrspeed);



  for(t=0;t<layers; t++)
    {

      for(j=0;j<50;j++) elstop[j][0]=elstop[j][1] = 0.0;
      switch (laycomp[t]) {
      case 2 :
	  if(nums[t][0]>nums[t][1]) {pie=nums[t][1];suu=nums[t][0];} else {pie=nums[t][0];suu=nums[t][1];}
	  nim=0;
	  for(n=1;n<=pie;n++) 
	      if(n>nim) nim=n;
	  printf("Two component layer, ratio=(%d:%d)\n",nums[t][0],nums[t][1]);
	  break ;
      case 1: 
	printf("One component layer\n");
	  break ;
      case 3 :
	printf("\n multi component layer\n");
	  break ;
      default :  /* pour les autres valeurs de i */
       break ;
     }

      if(laycomp[t]==1) {lukum=1;}
      else if (laycomp[t]==2){ lukum=2;}
      else {
	lukum=ntyp[t];
	printf("lukum %d \n",lukum);
      }


      for(n=0;n<lukum;n++) {
      sprintf(target, " %d ", multi_Z[t][n]); 
	printf("Target Z=%d\n",multi_Z[t][n]);
	if(laycomp[t]==2){
	  sprintf(zbl,"%s %s %s %s %s %s",inizbl,unit2,projectile,target,speeddat,end);
	}
	else{
	  sprintf(zbl,"%s %s %s %s %s",inizbl,projectile,target,speeddat,end);
	}
	printf("Executing ZBL command:\n  %s\n",zbl);

	system(zbl);
	gp = fopen("temp.in", "rt");
	for(j=0;j<50;j++) {
	  m=fscanf(gp, "%g %g", &speed, &force);
	  /*	  printf ("%g %g \n", speed, force);*/
	  if (m<2) {
	    printf("Failed reading elstop data from zbl output file.\n");
	    printf("Not creating elstop file. Check the zbl command.");
	    elstoperror=True;
	    break;
	  }
	  if(n<1){
	    elstop [j] [0]= elstop [j] [0] + (speed * 2187691.417);
	  }
	  if(laycomp[t]==1){ 
	    elstop [j] [1]= elstop [j] [1] + (force * 100.0); /*force is 1 keV/nm = 100 eV/Å -> eV/Å*/ 
	  }
	  else if(laycomp[t]==2) {
	    elstop [j] [1]= elstop [j] [1] + (nums[t][n]/nim)*force; /*sigma_ab=num(a)*sigma_a+num(b)*sigma_b*/
	  }  /*unit is eV/(1e15 atoms/cm²)*/
	  else {
	    rat=(float)nums[t][n]/numsa[t];
	    /*	    printf("%d %d %g \n",nums[t][n],numsa[t],rat);*/
	    elstop [j] [1]= elstop [j] [1] + force*100*rat; 
	    /*	    printf("elstop 0 1 %g %g \n",elstop [j] [0] ,elstop [j] [1] );*/
}
	}
	fclose(gp);
	system("rm temp.in"); 
      }
      if(laycomp[t]==2) {
	/* We have to apply Bragg's rule and therefore we must know the atomic density of the
	   compound */
	adens=1/(a*b*c); /*Å^-3*/
	adens*=(1e8)*(1e8)*(1e8);  /*cm^-3*/
	adens*=numsa[t];  /*atoms*cm^-3, numa=4 for fcc, so DO MORE CALC!!*/
	adens*=1e-23;  /*(1e23 atoms/cm^3)*/
	printf("What's the atomic density of this substance (1e23 atoms/cm^3) [%g]: ",adens);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &adens);
      }
      /*compound force=eV/(1e15 atoms/cm²), should be 1eV/Å
	eV/(1e15 atoms/cm²)*(1e23 atoms/cm^3)=1e8 eV/cm=eV/Å*/

      if (!elstoperror) {
	strcpy(fname,file2);
	if(multi) {sprintf(fname,"%s.%d",fname,t+1);}
	if((ep=fopen(fname, "rt"))!=NULL) {
	  printf("\n");
	  printf("The default elstop file (%s) already exists\n", fname); 
	  printf("Where do you wish to save the electronic stopping [%s]: ",fname); 
	  if(fgets(buf,79,stdin)!=NULL)
	    sscanf(buf,"%s", fname);
	  fclose(ep);
	}

	ep = fopen(fname, "wt");
	if(laycomp[t]=1){
	  for(h=0;h<50;h++) fprintf(ep, "%g %g\n", elstop[h][0], elstop[h][1]);
	}
	else{
	  for(h=0;h<50;h++) fprintf(ep, "%g %g\n", elstop[h][0], adens*elstop[h][1]);
	}
	fclose(ep);
      }

    } /*multilayer elstop-loop*/

  printf("\n\n");
  printf("#----------------------------------------------------------------#\n");
  printf("#			The Parameter File			 #\n");
  printf("#----------------------------------------------------------------#\n");
  OK[0]='A';
  if(expert==1) {
    printf("\nDo you want to make a parameter file (Y/N) [Y]: ");
    if(fgets(buf,79,stdin)!=NULL) sscanf(buf,"%c",OK);  
  }
  if(OK[0]!='N' && OK[0] !='n') {
    for(;;) {
      printf("\nThe repulsive potential:");
      printf("\nThe use of reppot-files is strongly recommended!");
      printf("\nIf you don't have a reppot-file, the ZBL-85 universal stopping power");
      printf("\nwill be used instead.");
      printf("\n\nDo you have a reppot-file for interactions between the projectile and");
      OK[0]='A';
      printf("\n");
      for(n=1;n<=donetypes;n++) {
	printf("    %s (Y/N) [N]: ", p_table[Z[n]-1].symbol);
	if(fgets(buf,79,stdin)!=NULL) sscanf(buf,"%c",OK);
	if(OK[0]=='Y' || OK[0]=='y')
	  reptyp[n]=2;
	else
	  reptyp[n]=szbl=1;
	OK[0]='A';
      }
      if(szbl==1 && expert==1)
	OK[0]='Y';

      if(OK[0]=='Y' || OK[0] =='y') {
	printf("\nZBL-85 parameters\n");
	printf("a\t b\n");
	for(k=0;k<6;k++) {
	  printf("%g\t %g\n", reppot[k][0], reppot[k][1]);
	}
	for(;;) {
	  OK[0]='A';
	  printf("Do you want to change the parameters (Y/N) [N]: ");
	  if(fgets(buf,79,stdin)!=NULL)
	    sscanf(buf,"%c",OK);
	  if(OK[0]=='Y' || OK[0]=='y') {
	    printf("\nPlease enter the a-parameters\n");
	    for(k=0;k<6;k++) {
	      printf("[%g] :", reppot[k][0]);
	      if(fgets(buf,79,stdin)!=NULL)
		sscanf(buf, "%g", &reppot[k][0]);
	
	    }
	    printf("Then the b-parameters\n");
	    for(k=0;k<6;k++) {
	      printf("[%g] :", reppot[k][1]);
	      if(fgets(buf,79,stdin)!=NULL)
		sscanf(buf, "%g", &reppot[k][1]);
        
	    }
	  }

	  else
	    break;
	}
      }
      if(multi!=1) {
	pakkot=0;
	OK[0]='A';
	printf("\n0.  Is the target polycrystalline (Y/N) [N]: ");
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%c",OK);
	if(OK[0]=='Y' || OK[0] =='y') {
	  help=50;
	  l=0;
	  Fiimax=360;
	  Thetamax=180;
	  for(;;) {
	    printf("\n    Please enter the average size of the domains (Å, ? for help): "); 
	    if(fgets(buf,79,stdin)!=NULL) {
	      if(buf[0]!=63) {
		sscanf(buf,"%g", &polysize);
		l=1;
	      }
	      else
		givehelp(help);
	    }
	    if(l==1) 
	      break;   
	  }
	  printf("\n    Please enter the standard deviation of domain size (Å): ");
	  scanf("%g", &polysizesigma);
	  printf("\n    Please enter the max. angle between domains (degrees, if 90 random orient.): ");
	  scanf("%g", &delta);
	  pakkot=1;
	  if(fgets(buf,79,stdin)!=NULL) {
	  }
	}
      }
		
      printf("\n1.  Please enter the number of projectiles calculated [%d]: ", Ncalc); 
      if(fgets(buf,79,stdin)!=NULL)
	sscanf(buf,"%d", &Ncalc);
      if(expert==1){
	printf("\n  Please enter the experimental Dose to compare with(realrange.out) [%g]: ", Dose);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &Dose);


	printf("Do you want to use sputtering modelling?\n");
	printf("(Y/N) [N]: ");
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%s",OK);
	if(OK[0]=='Y' || OK[0]=='y')
	  sputtering=1;

	printf("\n Please enter the amount of eroded surface at start [%g]: ", sput_surf);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &sput_surf);
      }

      if(pakkot==1) {
	Trange=1;
	rtype=2;
      }
      else {
	printf("\n2.  Please enter type of range\n");
	printf("    0=Z coordinate, 1=projection to the orig. direction of recoil atom\n"); 
	printf("    2=chord, 3=actual path of recoil atom [%d]: ", Trange); 
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%d", &Trange);
	if(Trange!=1) {
	  printf("\n3.  Please enter type of calculation\n");
	  printf("    1=normal, 2=negative ranges allowed [%d] : ", rtype);
	  if(fgets(buf,79,stdin)!=NULL)                                                             
	    sscanf(buf,"%d", &rtype);
	}
	else
	  rtype=2;
      }
      OK[0]='A';
      printf("\n4.  The angles are chosen randomly between the min. & the max.\n");
      printf("    Please enter the azimuth (Fii) angle range [%g - %g]:", Fii0, Fiimax); 
      if(fgets(buf,79,stdin)!=NULL) {
	n=sscanf(buf,"%g %c %g", &Fii0, OK, &Fiimax);
	if(n<3 && polysize==0) Fiimax=Fii0;
      }
      OK[0]='A';
      printf("    Please enter the polar (Theta) angle range [%g - %g]:", Theta0, Thetamax); 
      if(fgets(buf,79,stdin)!=NULL) {
	n=sscanf(buf,"%g %c %g", &Theta0, OK, &Thetamax);
	if(n<3 && polysize==0) Thetamax=Theta0;
      } 
      if(expert==1) {
	printf("\n5.  Starting pos. of the recoil atom is also chosen randomly. \n"); 
	printf("    X- and Y-positions are measured in unit cells and unless you want to\n"); 
	printf("    emphasise some region, you should choose the values so that max = min +1\n"); 
	printf("    Please enter the x-range [%g - %g]:", minx, maxx);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g - %g", &minx, &maxx);
	printf("    Please enter the y-range [%g - %g]:", miny, maxy);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g - %g", &miny, &maxy);
	printf("    Please note that the Z-position must be small enough that \n");
	printf("    the recoil atom is not created inside the lattice!\n");
	printf("    Please enter the z-range (Å) [%g - %g]:", minz, maxz);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g - %g", &minz, &maxz);
      }
      OK[0]='N';
      if(expert==1) {
	printf("\nDo you wish to change any of the parameters 0-5 (Y/N) [N]: ");
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf," %c", OK);
      }
      if(OK[0]!='Y' && OK[0]!='y')
	break;
    }
    for(;;) {
      if(expert==1) {
	printf("\n6.  Moving limit: the fraction of the unit cell that must\n");
	printf("    always remain in front of the recoil atom; if the atom moves\n");
	printf("    past this limit, fresh crystal shall be created\n");
	printf("    ! The moving limits should be larger than the potential cut-off !\n");
	printf("    ! Do not change them unless you're sure you know what you are doing !\n");
	printf("    Please enter the x-limit [%g]: ", lim[0]);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &lim[0]);
	printf("    Please enter the y-limit [%g]: ", lim[1]);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &lim[1]);
	printf("    Please enter the z-limit [%g]: ", lim[2]);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &lim[2]);
      }
      printf("\n7.  Please enter the temperature of the initial state(K) [%g]: ", Tini); 
      if(fgets(buf,79,stdin)!=NULL)
	sscanf(buf,"%g", &Tini);
      printf("\n");
      printf("8.  The initial state calculation: \n");
      printf("    Please note that the Debye temp. depends on the real sample temperature\n"); 
      for(n=0;n<layers;n++) {
	printf("\n");
	printf("    Please enter the Debye temperature for %d. layer\n", n+1); 
	printf("    (K, 0=use an attractive potential instead) [%g]: ", dtemp[n]); 
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf, "%lf", &dtemp[n]);
	if(dtemp[n]<=0)  {
	  printf("    Which attractive potential should be used?\n");
	  printf("    1 Mazzone 2 Morse : ");
	  scanf("%d", &atrpot);
	  if(fgets(buf,79,stdin)!=NULL) {
	  }
	}   
      }
      if(expert==1) {
	printf("\n9.  How many iterations to determine the initial configuration [%d]: ",Ninistate);
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%d", &Ninistate);
      }
      l=0;
      help=48;
      if(dtemp[n]>0) {
	for(;;) {
	  printf("\n    Please enter the time of the initial state calculation\n");
	  printf("    (fs, ? for help) [%g]: ", Timeini); 
	  if(fgets(buf,79,stdin)!=NULL) {
	    if(buf[0]!=63) {
	      sscanf(buf,"%g", &Timeini);
	      l=1;
	    }
	    else
	      givehelp(help);
	  }
	  if(l==1)
	    break;
	}
      }	
      seed=getpid(); /* the seed is created from pid & clock */
      if(seed<=500)  
	seed=seed*seed;
      seed=walltime%seed;
      if(expert==1) {
	printf("\n10. Please enter a seed for the random number generator");
	printf("\n    (integer) [use pid & clock (%d)]: ", seed); 
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%d", &seed);
      }
      for(;;) {
	if(dtemp[n]>0 && expert==1) {
	  printf("\n    (NEIGBOUR LIST CUT-OFF > POTENTIAL CUT-OFF)");
	  printf("\n11. Please enter the NEIGHBOUR LIST construction cut-off radius\n"); 
	  printf("    to be used in INICALC.C(Å) [%g]: ", Rmini); 
	  if(fgets(buf,79,stdin)!=NULL)
	    sscanf(buf,"%g", &Rmini);
	  printf("    Please enter the POTENTIAL cut-off radius\n");
	  printf("    to be used in INICALC.C(Å) [%g]: ", R0ini);
	  if(fgets(buf,79,stdin)!=NULL)
	    sscanf(buf,"%g", &R0ini);
	  printf("\n");
	}
	if(Rmini >=  R0ini)
	  break;
      }
      if(expert==1) {
	help=49;
	l=0;
	for(;;) {
	  printf("\n15. By what coefficient do you want to scale the elstop (? for help) [%g]: ", elscale);
	  if(fgets(buf,79,stdin)!=NULL)  {
	    if(buf[0]!=63) {
	      sscanf(buf,"%g", &elscale);
	      l=1;
	    }   
	    else
	      givehelp(help);
	  }
	  if(l==1)
	    break;
	}
	printf("    By what coefficient do you want to scale the repulsive pot. [%g]: ", repscale); 
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf,"%g", &repscale);
      }
      OK[0]='N';
      if(expert==1) {
	printf("\nDo you want to change any of the parameters 6-15 (Y/N) [N]: ");
	if(fgets(buf,79,stdin)!=NULL)
	  sscanf(buf," %c", OK);
      }
      if(OK[0]!='y' && OK[0]!='Y')
	break;
    }


    if((hp=fopen(file3,"rt"))!=NULL) {
      printf("\n");
      printf("The default parameter file (%s) already exists\n",file3); 
      printf("Where do you want to save the parameters[%s]: ", file3); 
      if(fgets(buf,79,stdin)!=NULL)
	sscanf(buf,"%s", file3);
      fclose(hp);
    }
    hp=fopen(file3,"wt"); 
    for(n=0;n<=donetypes;n++) {
      fprintf(hp,"type[%d].Z:= %d\n", n, Z[n]);
      fprintf(hp,"type[%d].m:= %g\n", n, mass[n]);
    }
    fprintf(hp,"physical->Dose:= %g\t# Physical Dose for experimental comparison\n",Dose);	
    fprintf(hp,"physical->Damage:= %d\t# Flag for Damage modelling (0=no,1=yes)\n",Damage);
    fprintf(hp,"physical->scaling:= %d\t# Flag for dimensional scaling during the simulation\n",scaling);
    fprintf(hp,"			      according to the volume in the coords.in.# files\n");
    fprintf(hp,"physical->bscaling:= %d\t# Flag for dimensional scaling at the start of simulation\n",bscaling);
    fprintf(hp,"			      according x_i/box.movelim.i ratio\n");
    fprintf(hp,"physical->sputtering:= %d\t# Flag for sputtering modelling\n",sputtering);
    fprintf(hp,"physical->sput_surf:= %g\t# amount of eroded surface in Ångströms\n",sput_surf);

    fprintf(hp, "reccalc->E0:= %g\t# Recoil energy \n", E0);
    fprintf(hp, "reccalc->Atype:= 0\t# Recoil atom type\n");
    fprintf(hp,"reccalc->Ncalc:= %d \t# Number of histories\n", Ncalc);
    fprintf(hp, "reccalc->Trange:= %d\t# Range type: 0=sz, 1=proj, 2=chord, 3=path\n", Trange); 
    fprintf(hp, "reccalc->type:= %d\t# Recoil calculation type: 1=normal,\n", rtype); 
    fprintf(hp,"                          2=negative range allowed\n");
    fprintf(hp, "reccalc->Fii0:= %g\t# Azimuth angle of the recoil atom\n", Fii0);
    fprintf(hp, "reccalc->Fiimax:= %g\t# is chosen randomly from [Fii0,Fiimax]\n", Fiimax); 
    fprintf(hp, "reccalc->Theta0:= %g\t# Polar angle of the recoil atom\n", Theta0); 
    fprintf(hp, "reccalc->Thetamax:= %g\t# is chosen randomly from [Theta0,Thetamax]\n\n\n",Thetamax); 
    if(pakkot==1) {
      fprintf(hp, "reccalc->Polysize:= %g \t# Polycrystalline settings\n",polysize);
      fprintf(hp, "reccalc->Polysizesigma:= %g\n", polysizesigma);
      fprintf(hp, "reccalc->Polydeltaangle:= %g\n", delta);
    }
    fprintf(hp, "reccalc->Startmin.x:= %g\t# The starting pos. of the recoil\n", minx*multi_a[0]);
    fprintf(hp, "reccalc->Startmax.x:= %g\n", maxx*multi_a[0]);
    fprintf(hp, "reccalc->Startmin.y:= %g\n", miny*multi_b[0]);
    fprintf(hp, "reccalc->Startmax.y:= %g\n", maxy*multi_b[0]);
    fprintf(hp, "reccalc->Startmin.z:= %g\n", minz);
    fprintf(hp, "reccalc->Startmax.z:= %g\n\n", maxz);
    if(atrpot!=0)
      fprintf(hp, "pot->attr.type:= %d\n", atrpot); 
    for(n=1;n<=donetypes;n++)
      fprintf(hp, "pot->rep.type[0][%d]:= %d\n\n",n, reptyp[n]);
    if(szbl==1) {
      for(n=0;n<6;n++) {
	fprintf(hp, "pot->rep.ai[%d]:= %g\n", n, reppot[n][0]);
	fprintf(hp, "pot->rep.bi[%d]:= %g\n", n, reppot[n][1]);
      }
    }

    fprintf(hp, "layers->Nlayers:=%d\n", layers);
    fprintf(hp, "layers->Natomsmax:= %d\n", 20*lkm);
    fprintf(hp, "layers->Ninistate:= %d\n\n", Ninistate);
    if(multi==1)  {
      for(t=0;t<layers;t++) {
	fprintf(hp, "layers->layer[%d].Probmode:= 0\n", t);
	fprintf(hp, "layers->layer[%d].minz:= %g\n", t, laydim[0][t]);
	fprintf(hp, "layers->layer[%d].maxz:= %g\n", t, laydim[1][t]);
	if(layers>1)
	  fprintf(hp, "layers->layer[%d].ltype:= %d\n", t, t);
	fprintf(hp, "layers->layer[%d].debyet:= %g\n\n", t, dtemp[t]);
      }
    }
    else {
      fprintf(hp, "layers->layer[%d].debyet:= %g\n\n", 0, dtemp[0]);
	    
    }
    fprintf(hp, "box->Max.x:= %g\t#The size of the simulation cell\n", boxx);
    fprintf(hp, "box->Max.y:= %g\n", boxy);
    fprintf(hp, "box->Max.z:= %g\n\n", boxz);
    lim[0]=lim[0]*multi_a[0];
    lim[1]=lim[1]*multi_b[0];
    lim[2]=lim[2]*multi_c[0];
    fprintf(hp, "box->movelim.x:= %g # Shift used in generating fresh crystal\n", lim[0]); 
    fprintf(hp, "box->movelim.y:= %g # in front of the recoil.\n", lim[1]);
    fprintf(hp, "box->movelim.z:= %g\n\n", lim[2]); 
    fprintf(hp, "physical->Tini:= %g # Temperature of the inistate calculation\n", Tini);
    if(dtemp[0]<=0)
      fprintf(hp, "time->Timeini:= %g\n", Timeini);
    fprintf(hp, "gen->seed:= %d\n\n", seed);
    if(dtemp[0]<=0) {
      fprintf(hp, "potcrit->Rmini:= %g # cut-off radiuses for inistate calc.\n", Rmini);
      fprintf(hp, "potcrit->R0ini:= %g\n", R0ini);
    }
    fprintf(hp, "potcrit->Rmrec:= %g # cut-off radiuses for recoil calc.\n", Rmrec);
    fprintf(hp, "potcrit->R0rec:= %g\n\n", R0rec);
    fprintf(hp, "elstop->scale:= %g # Electronic stopping power scaled by this\n", elscale);
    fprintf(hp, "pot->rep.scale:= %g # Repulsive potential scaled by this\n", repscale);
    fclose(hp);
  }
  printf("\n");
  printf("#----------------------------------------------------------------#\n");
  printf("#			The MDH Script				 #\n");
  printf("#----------------------------------------------------------------#\n");
  OK[0]='A';
  printf("\n\nDo you want to make an MDH-script (Y/N) [N]: ");
  if(fgets(buf,79,stdin)!=NULL)
    sscanf(buf,"%c",OK);
  if(OK[0]=='Y' || OK[0] =='y') {
    OK[0]='A';
    printf("\nShould MDH calculate the deposited energies tables [N]: ");
    if(fgets(buf,79,stdin)!=NULL)
      scanf(buf, "%c",OK);
    if(OK[0]=='Y' || OK[0] =='y')
      strcat(mdh, " -FD");
    OK[0]='A';
    printf("Should MDH calculate the stopping tables [N]: ");
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf,"%c",OK);
    if(OK[0]=='Y' || OK[0] =='y')
      strcat(mdh, " -S");
    OK[0]='A';
    printf("Should MDH calculate the average velocities tables [N]: ");   
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf, "%c",OK);
    if(OK[0]=='Y' || OK[0] =='y')
      strcat(mdh, " -AV");
    OK[0]='A';
    printf("Should MDH do statistics on E(z) [N]: ");
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf,"%c",OK);
    if(OK[0]=='Y' || OK[0] =='y')
      strcat(mdh, " -Es");
    help=60;
    l=0;   
    for(;;) {
      OK[0]='A';
      printf("Do you want debug output (? for help) [N]: ");
      if(fgets(buf,79,stdin)!=NULL) {
	if(buf[0]!=63) {
	  sscanf(buf,"%c",OK);
	  l=1;
	}
	else
	  givehelp(help);
      }
      if(l==1)
	break;
    }
    if(OK[0]=='Y' || OK[0] =='y') {
      strcat(mdh," -d");
      OK[0]='A';
      printf("Do you want maximum debug [N]: ");
      if(fgets(buf,79,stdin)!=NULL)
	sscanf(buf,"%c",OK);
      if(OK[0]=='Y' || OK[0] =='y') {
	strcat(mdh,"m");
      }
    }
    printf("Which param-file should be used [%s] : ", file3);
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf,"%s",file3);
    strcat(mdh," -pa ");
    strcat(mdh,file3);
    if(szbl!=1 && comp!=2 && multi!=1) {
      printf("Which reppot-file should be used [%s] :", file4);
      if(fgets(buf,79,stdin)!=NULL) {
	sscanf(buf,"%s",file4);
	strcat(mdh," -rp ");
	strcat(mdh,file4);
      }
    }
    printf("\n");
    printf("Nice is a UNIX command that prevents your program from hogging all\n");
    printf("the resources and thus the other users will still be NICE to you\n");
    printf("How nice [8]: ");
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf,"%d", &nice);
    sprintf(nicei, "%d ", nice);
    strcat(nices, nicei);
    strcat(nices, mdh);
    printf("Do you want to run MDH in the background (Y/N) [Y]: ");
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf, " %c", OK);
    if(OK[0]!='N' && OK[0] !='n')
      strcat(nices," &");
    printf("Where do you want to save this script [%s]: ", file5);
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf,"%s", file5);
    fp=fopen(file5, "wt");
    fprintf(fp,"%s", nices);
    fclose(fp);
    OK[0]='A';
    strcat(chmod,file5);
    system(chmod);
    if(szbl!=1 && multi!=1)
      printf("\nTo run mdh, you must have a reppot-file available!\n");
    printf("Do you wish to run the script now [N] : ");
    if(fgets(buf,79,stdin)!=NULL)
      sscanf(buf, " %c", OK);
    if(OK[0]=='Y' || OK[0] =='y')
      system(nices);
  }
  else {
    printf("\nYou can now start the range simulation with the command mdh\n");
  }
  printf("\n Remember that you are now using ZBL-el.stopping model.\n");
  printf("If you are using Si as a lattice material and you want\n");
  printf("to use PENR/BK-models, look at the mdh3_* - directories for examples!\n");
  printf("\nMay the force be with you!\n");
}


void givehelp(int help) {
  char ch;
  FILE *fp;
  char endhelp[]="help.txt";
  char helps[70];
  strcat(helps,DATADIR);
  strcat(helps,endhelp);

  fp = fopen(helps, "rt");
  for(;;) {
    while((ch=fgetc(fp)) != 36) {
    }
    ch=fgetc(fp);
    if (ch == help) {
      break;
    }
  }
  while((ch=fgetc(fp))!= 35) {
    printf("%c", ch);
  }
  printf("\n");
  fclose(fp);
}

