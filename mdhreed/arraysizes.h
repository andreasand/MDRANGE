/* Some of these could later be changed to use malloc ! */


/* maximum number of possible atom types   */
/* Must be greater than 1, for some reason */

#define MAXTYPE 10

/* maximum number of allowed layers */
#define ORIGMAXLAYER 4
EXTERN int MAXLAYER;

#define MAXSPLITBORDERS 6
#define MAXSPLITATOMS 64  /*2^MAXSPLITBORDERS*/

/* maximum number of allowed different inistate calcs for one layer */
#define ORIGMAXINISTATE 10
EXTERN int MAXINISTATE;

/* Maximum number of atoms allowed */
#define ORIGMAXATOM 1000 
EXTERN int MAXATOM;

/* Maximum number of amorphization states allowed */
#define ORIGMAXAMORPH 20 
EXTERN int MAXAMORPH;

#define MAXN 16 /*maximum number of firsov partners*/

/* Maximum number of unitboxes in z-direction in damage calculation JARKKO*/
#define ZBOXES 10000

/* Maximum numbers of neighbours possible */
#define NGBRSTEP 10 
EXTERN int MAXNGBR;

/* Max number of exponential terms in potential */
#define MAXPOTEXP 6  

#define MAXRHO 128  /*size of the charge density table, if needed*/
#define MAXELSTOP 1000  /*max number of zbl points*/
#define MAXSTOPP  200   /*max number of cdensity and stopping points*/
#define MAXRECOIL 100000

#define DEPENSIZE 5000
#define STOPSIZE 100
#define AVVELSIZE 250

#define RECSPECSIZE 400

/* Array sizes for statistics of E(z), fii(z) and theta(z) */
/* fii varies between 0 and 360, theta 0 and 180 */
#define ESTATZSIZE 50
#define ESTATESIZE 50
#define ESTATFIISIZE 72 
#define ESTATTHETASIZE 36

#define MAXMORSE 20
#define RARRAYSIZE 4000
#define VELSTORESIZE 15

/* Size of arrays for read in repulsive potential data */

#define REPPOTMAX 502

#define PI 3.1415927

#define a_b 0.529177343575629789 /*Bohr radius for converting to atomic units*/


