/* 
  All structures used in many places in the program and passed from 
  subroutine to other. References are to Kai Nordlunds gradu 
*/



/* 
  UNIT structures
*/
struct unit {
   real energy,length,force;
   real tempconv;
   real mass[MAXTYPE];
   real time[MAXTYPE],vel[MAXTYPE],acc[MAXTYPE];
};

/* 
  Miscellaneous SMALL structures
*/
struct vvector { real x,y,z; };
typedef struct vvector vector;
typedef struct vvector point;
struct polar  { real r,theta,fii; } polar;

struct gen {
   int seed;    /* Seed for random number generator */
};

/* 
   Structures that deal with ATOMs
*/
struct atom {   /* basic atom properties (coordinates mainly) */
   real x,y,z;
   real vx,vy,vz;
   real ax,ay,az,a1x,a1y,a1z,a2x,a2y,a2z;    /* accelerations for different times */
   int type;    /* The type structure holds the other physical properties */
};
struct atomex {  /* Miscellaneous properties */
   real pot;   /* Total potential the atom experiences */
   real force; /* Sum of individual force magnitudes */
};
struct atomflag {   /* flags etc. */
   logical moving;
   logical recoil;
   logical energetic;   /* Has the atom received energy from recoil ? */
   logical interact;    /* Does the atom experience an interaction this time step */
};
struct type {
   int Z;     /* Atomic number */
   real m;    /* REAL mass in u. In internal units, masses are all 1 */
   real q;    /* charge */
   real prob; /* Probability of occuring, if 1.0 no effect */
};
struct time { /* Because of the unit system (mass=1) these depend on atom type */
   real dtreal,dtmaxreal,prevdtreal; 
   real dt[MAXTYPE],prevdt[MAXTYPE],dtratio[MAXTYPE];   /* Internal time steps */
   real tottime;    /* Total real time */
   real kt,Kt;      /* criterion for time step modification */
   real Et;         /* criterion for time step modification during strong collisions */
   real Dt,Dtini;          /* Original time step in recoil and inistate calc. */
   real Dtmax;             /* Maximum time step for tstep increasing */
   real Timeini;      /* Length of inistate calc. in fs */
};
struct iniat {   /* Ignore accelerations in storing initial atom coords */
   real x,y,z,vx,vy,vz; 
   int type;
};

/* 
   Structures that deal with the structure of the sample
*/
struct layer {
  int fulldensity; /*flag for using full charge density table*/
   int elstop;    /*flag for elstop mode to use in this layer: 0=zbl, 1=puska for ion-Si,2=puska for ion-Si+ion-ion JARKKO*/
   real minz,maxz;     /* Depth limits at which layer may be on */
   real density;       /* layer density */
  real surfbind;       /*binding energy for the surface atoms in this layer*/
   real SIdensity;     /* layer density in SI units*/
   real ndensity;      /* number density, i.e. nat / volume of cell */
   real Zmean;         /* Weighted mean Z of substrate atoms */
   int ltype;          /* Type is number of read-in coordinate file coords.in */
  int natoms[ORIGMAXAMORPH];   /* Number of atoms in layer in each am. state*/
  point boxsize[ORIGMAXAMORPH];  /*boxsizes for each am. state+layer*/
  point unitboxsize[ORIGMAXAMORPH];  /*real unitbox sizes for the layer*/
  point front[ORIGMAXAMORPH];  /*front coords for layer change;*/
  real Dose[ORIGMAXAMORPH];
   real Timeinior;     /* Override for time of inistate calc. for this layer */
   real Probability;   /* Probability that this layer should be selected at 
			  the correct depth, between 0 and 1, 1.0 default */
   int Probmode;       /* Mode of selection: 0: select at every new moveslice call */
                       /* 1: select only first time minz is crossed, 2: as 0, but */
		       /* use Thickness parameter */
   real Thickness;     /* Thickness of the layer once selected, if Probmode = 2 */
   logical crossed;   /* Has the layer minz limit been crossed during this recoil ? */
   logical selected;   /* Has the layer been selected during this recoil ? */
   double debyet;      /* Debye temperature, 0 in N/A */
  int Namorph;         /*number of am. states */
  real steadydens;     /*steady electron density for spherical symmetry elstop models*/
  int shiftmode;       /* Should cell be shifted over periodic boundaries */
};
struct layers {
   int Nlayers;
   int Ninistate;    /* how many inistate calcs for each layer ? */ 
   int Natomsmax;    /* Max number of atoms for array allocation */
   logical onelayer;
   struct layer *layer;
};


/* 
  Structures that deal with POTENTIALS and STOPPING
*/
struct potcrit {    /* cut-off radiuses and other criteria */
   real r0,rm; /* rm is outer radius, r0 the true potential cut-off radius */ 
   real R0rec,Rmrec,R0ini,Rmini; /* recoil and inistate calc values */
   int fatoms;  /* No. of atoms closer than Firsov cut-off */
   real frmin[MAXN]; /* Dist. to nearest Firsov-partners */  
   int fin[MAXN]; /* Indices of Firsov-partners */
   real rfirsov; /* Cut-off for the Firsov model */
};
struct reppot {  /* Repulsive potential vars */
   int type[MAXTYPE][MAXTYPE];       /* Potential type: 1 ZBL-type, 2 spline interpolate file */
   double k_e,fact[MAXTYPE][MAXTYPE];     /* Electric constant and k_e*Z_1*Z_2 */
   double au[MAXTYPE][MAXTYPE];           /* a factor, eq. (2.29) page 19 */
   double ai[MAXPOTEXP],bi[MAXPOTEXP];    /* exponential terms for readin */
   double a[MAXTYPE][MAXTYPE][MAXPOTEXP],b[MAXTYPE][MAXTYPE][MAXPOTEXP];  /* exponential terms */
   real scale;                            /* Rep. potential scaling factor */
   int nterms;
};
struct attrpot {  /* Attractive potential factors, whatever */
   int type;    /* 1 for Mazzone, 2 for Morse */
   real ljeps[MAXTYPE][MAXTYPE],ljsigma[MAXTYPE][MAXTYPE];

   real mazD,mazalpha,mazr1,mazd,mazK,mazr2;
   
   int morsetype;   /* See potential.c */
   real morsealpha[MAXMORSE],morser1[MAXMORSE],morseD[MAXMORSE];

};
struct pot {
   struct reppot rep;
   struct attrpot attr;
};
struct elstop {
  int Firsov; /*flag*/
  real totalfirsov; /*total force due to the firsov model*/
   int n[ORIGMAXLAYER];       /* Number of readin points for ZBL file for each layer*/
   int nc[MAXTYPE];           /* Number of readin points for charge density file for each atom in layers*/
   int ns;                    /* Number of readin points for stopping_density file for ion*/
   int mode;                 /* mode: if 1 include straggling of el. stop. */
   real scale;               /* Elstop scaling factor at E=0*/
   real dscale;              /* Linear modification */
   real sqscale;             /* Quadratic mod., final scale=scale+dscale*E+sqscale*E*E */
   real v[ORIGMAXLAYER][MAXELSTOP],Setable[ORIGMAXLAYER][MAXELSTOP]; /*Number of current layer should be added*/
  real cdist[MAXTYPE][MAXSTOPP],celdens[MAXTYPE][MAXSTOPP];  /*el.density vs. r*/
  real strs[MAXSTOPP],stopp[MAXSTOPP];     /*Se vs. rs (for puska) or G*Slin vs. rs (for BK)*/
  real slim,iperslim;   
   real Se;
  real *rho;  /*if full charge density file provided, it will be loaded into giant table, rho*/
  real dx,dy,dz;  /*distance from point to point in previous 3d-box (size/128) */
  real uy,ux,uz;   /*ch.density box side lengths*/
  int correction;
  real dens[DEPENSIZE];
  real Sespc[DEPENSIZE];
  int penr; /*flag*/
};

/*
   Structures that deal with the RECOIL and the type of the recoil calc.
*/
struct recoil {   /* Data of an individual recoil during simulation */
   int n;        /* Number of current recoil calculation */
   int number;  /* number of the recoil particle (atom index) */
   int attype;    /* It's atom type */
   real energy,vel,acc;  /* calculated for every time step */
   real x0,y0,z0;           /* It's original (starting) position in the current domain */
   real vx0,vy0,vz0,vel0;   /* It's original velocity vector in the current domain*/
   real x,y,z;            /* It's position in the simulation cell */
   real vx,vy,vz;         /* It's current velocity vector */
   real xp,yp,zp;         /* It's previous position in the simulation cell */
   real eprev,vprev;      /* It's previous energy and velocity*/
   real dx,dy,dz,dr,drp;  /* last motion */
   real sx,sy,sz;      /* Total position */
  
   real px,py,pz;      /* sx + sum over polycryst. domains */
   real r_path;        /* range along the path (sum of dr:s) */
   real r_chord;       /* chord range */
   real r_proj;        /* Range projected to original velocity vector */
   real range;	       /* Current range selected by reccalc->Trange */
   real theta,fii;     /* Original theta, fii */
   logical beeninside;     /* Has the recoil been inside (z > box->movelim.z) yet ? */
   logical realbeeninside;     /* Is the recoil inside the real surface ? */
   logical turnedaround;   /* Has the z velocity been negative yet ? */
   logical incollision;    /* Is the recoil in a violent collision right now ? */
  logical touched;
  real ddx,ddy,ddz; /*total error from boxscaling (is negligible)*/
};

struct maxima {
  real vel,acc;
  int vel_attype,acc_attype;

};

struct reccalc {   /* Data of the whole recoil calculation type etc. */
   int type;       /* Recoil calc. type, see recoil.c */
   int Atype;      /* Recoil atom type */
   int Ncalc;      /* How many recoils should be calculated ? */
   int Trange;     /* Which range should we use: 0 sz, 1 proj, 2 chord, 3 path */
   logical Estat;  /* Should we do statistics of recoil E, theta and fii (z) ? */
   real Estatzmax; /* Which z should we use as max in a Estat calc. ? */
   real E0;        /* Initial recoil energy */
   real v0;        /* Initial velocity, calculated from E0 */
   real Emin;      /* Minimum energy that the particle has */
   real Zmax;      /* End calculation at z=zmax */
   point Startmin,Startmax;
   real Theta0,Fii0,Thetamax,Fiimax;
   real Polysize,Polysizesigma;  /* Assume gaussian distributed grain size */
   real Polydeltaangle;          /* Maximum change in direction, in degrees */
   logical poly;                 /* True if Polysize != 0.0 */
   int Forceseed;                /* Force this seed for recoil calc. */
   int Outputstep;               /* Interval for coordinate output */
   logical change; /* If any type[].prob != 1.0, change set to true */
   real recspeczmin,recspeczmax;  /* Z range for recspec calculation */
   real ETmax,ETmin;              /* Energy-range of depen angle analysis */
   real collisionsout;            /* Output collisions.txt file of all recscpec recoils? */
   int printz; 
   real stopDR; 
   int allowhugeenergy;
};

struct recresult {   
  double stopE,stopN,stopv,stopR,TstopE[STOPSIZE],TstopN[STOPSIZE];
  int istop,iTstop[STOPSIZE];

   real ranges[MAXRECOIL];   /* Table of final ranges */
   real range,straggling;    /* Mean range and straggling */
   double depennucl[DEPENSIZE],depenel[DEPENSIZE]; /* Deposited energy tables */

  double depennucl_box[ZBOXES],depennucl_boxtemp[ZBOXES];
  double surfFDn;   /*deposited energy in the surface layer for sputtering*/

   double sumFD,sumFDn,sumFDe,sumFDf;   /* Integrated value */
   double sumFDsput,sumFDtrans;  /* Energy deposited to sputtered 
				    and transmitted atoms */
   real dependz;                 /* Step used in depen tables */

   double nuclstop[STOPSIZE],elstop[STOPSIZE]; /* Nuclear and electronic stopping tables */
  double pnuclstop[STOPSIZE],pcurnuclstop[STOPSIZE]; /*JARKKO,for positive values*/
  int pinuclstop[STOPSIZE],picurnuclstop[STOPSIZE];  /*JARKKO,for positive values*/
   int inuclstop[STOPSIZE],ielstop[STOPSIZE]; 
   double curnuclstop[STOPSIZE],curelstop[STOPSIZE]; 
   int icurnuclstop[STOPSIZE],icurelstop[STOPSIZE]; 
   real stopdv,stopvmax;         /* Step and max vel used in stopping tables */ 

   double recspec[MAXTYPE][RECSPECSIZE];  /* Recoil spectrum arrays */
   int irecspec[MAXTYPE][RECSPECSIZE],irectot[MAXTYPE],irec0[MAXTYPE];
   double drecspec,rec0[MAXTYPE],rectot[MAXTYPE];
  int irecspec_Theta[MAXTYPE][RECSPECSIZE]; 

   double avvel[AVVELSIZE];      /* Average velocity (t) tables */
   int iavvel[AVVELSIZE];
   double curavvel[AVVELSIZE];   /* Average velocity (t) tables for current event */
   int icuravvel[AVVELSIZE];
   real avveltmax,avveldt;       /* Maximum time used in avvel tables */
   int nvalidavvel;  /* Number of cases that gave valid avvel calculation */

   int nrange;                   /* Number of recoils that gave a valid range */
   int stopway[8];               /* Table of number of recoils stopped different ways, see recoil.c */
};


struct poly { /* Internal data of polycrystallinity calc; readin data got from struct reccalc! */
   int ndomain;      /* Number of current domain */
   int ndomainstot;  /* Total number of domains for all recoils */
   real cursize;     /* Current grain size */
   real prevvel;     /* Recoil velocity in previous domain before boundary reached */
   vector prevv;     /* Recoil velocity vector -------------"------------------    */
   real range;       /* Polycrystalline range over previous domains */
   point recpos;     /* Position of recoil over previous domains */
   logical boundreached;    /* Did the recoil end because a boudary was reached ? */
};

/* 
  Structures dealing with the whole simulation cell
*/
struct box {        /* Simulation cell sizes etc. */
   point Min,Max,size;
   point movement;          /* Movement in recoil calc. */
   point movelim;           /* Limits for moving */
   logical periodic;
   logical warnoutside;
   int nomoveslice;      /* Allows preventing moveslize completely(1) or in z(2) */
   int Iniscaleint;      /* Time step interval between cm scalings in inistate calc. */
};
struct physical {    /* Physical constants of the simulation cell */
   real t,ekin,epot,etot;
   real p;
   real Tini,Trec;    /* Required temperatures */
   real maxv;
  
  int spreadFD;  /*flags ->*/
  int Damage;   
  int bscaling;  
  int scaling;   
  real Dose;   /*real dose*/
  real Dosescaling;  /*scaling factor for dose*/
  int sputtering;    /*flag*/
  real sput_surf;  /*amount of surface eroded, in Ongs.*/

  real deltasput_surf;  
  real Y; /*fixed value for yield if needed*/
  real sputFD; /*fixed value for surface depon. energy*/
  real sputdeltaV; /*value for the increased volume/original volume, because of the foreign ion in material usually [0-2]*/
  real rmin;  /* Minimum separation from recoil to atoms */
};

/*
  Structures containing I/O flags, names etc.
*/
struct file {
   char coord[120],param[120],elstop[120];
   char attrpot[120],reppot[120];
   FILE *fpco,*fppa,*fpel,*fpat,*fpre,*fpcdens,*fpstopp;
   logical readattrpot;

/* Output files. fpo=stdout and fpdb=stderr assumed */
   char outf[120],startf[120],rangef[120],depenf[120],stopf[120];
   char avvelf[120],debugf[120],recspecf[120];
  char recspec_Thetaf[120]; 
  char cdens[120],stopp[120],crho[120];
  FILE *fpo,*fps,*fpra,*fpde,*fpst,*fpav,*fpdb,*fprs,*fprst,*fprho,*fpcollisions;

   char estatef[120],estatff[120],estattf[120];  /* Estat files */
   FILE *fpee,*fpef,*fpet;

   int output;
  int movie; 
   logical fullcalc,printstart,printcoords,printcoordsi,printrec;
   logical printrange,printdepen,printstop,printavvel;
   logical recoilspec,printrecend,printdens;

  int layeratoms[MAXTYPE][ORIGMAXLAYER]; 
};

