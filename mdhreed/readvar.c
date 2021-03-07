#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* INTEGER READS */

int readi(varname,fp)
char *varname;
FILE *fp;
{
   int i=1;
   double d=0.0;
   int in;
   char b[2];

   readvar(varname,&in,&d,b,i,fp);
   return(i);
}
int readic(varname,v,v0,mess,fp,fps)
char *varname;
FILE *fp,*fps;
int *v,v0;
char *mess;
{
   int i=1;
   static int n;
   char b[2];
   double d;

   n=1;
   if (readvar(varname,v,&d,b,i,fp) >= 0) {
      if (mess[0]!='\0') fprintf(fps,"%s %s %d (read in)\n",mess,varname,*v);
   }
   else {
      *v=v0; n=0; 
      if (mess[0]!='\0') fprintf(fps,"%s %s %d (default)\n",mess,varname,*v);
   }
   return(n);
}

/* FLOAT (REAL) READS */

float readr(varname,fp)
char *varname;
FILE *fp;
{
   static int i,j;
   char b[2];
   double d;

   i=0;
   readvar(varname,&j,&d,b,i,fp);
   return((float) d);
}
int readrc(varname,v,v0,mess,fp,fps)
char *varname;
FILE *fp,*fps;
float *v,v0;
char *mess;
{
   int i=0,s,j;
   char b[2];
   double d;
   int readvar();
   static int n;

   n=1;
   s=readvar(varname,&j,&d,b,i,fp);
   *v=(float) d;

   if (s >= 0) {
      if (mess[0]!='\0') fprintf(fps,"%s %s %g (readin)\n",mess,varname,*v);
   }
   else {
      *v=(float) v0; n=0; 
      if (mess[0]!='\0') fprintf(fps,"%s %s %g (default)\n",mess,varname,*v);
   }
   return(n);
}

/* DOUBLE READS */

double readd(varname,fp)
char *varname;
FILE *fp;
{
   int i=0,j;
   char b[2];
   double d;

   readvar(varname,&j,&d,b,i,fp);
   return((float) d);
}
int readdc(varname,v,v0,mess,fp,fps)
char *varname;
FILE *fp,*fps;
double *v,v0;
char *mess;
{
   int i=0,n=1,s,j;
   char b[2];
   int readvar();

   s=readvar(varname,&j,v,b,i,fp);
   if (s >= 0) {
      if (mess[0]!='\0') fprintf(fps,"%s %s %.10lg (readin)\n",mess,varname,*v);
   }
   else {
      *v=v0; n=0; 
      if (mess[0]!='\0') fprintf(fps,"%s %s %.10lg (default)\n",mess,varname,*v);
   }
   return(n);
}

/* STRING READS */


char *reads(varname,fp)
char *varname;
FILE *fp;
{
   static char *s;
   int i=2;
   int j;
   double d;
   s=(char *) malloc(80);
   readvar(varname,&j,&d,s,i,fp);
   return(s);
}
int readsc(varname,v,v0,mess,fp,fps)
char *varname;
FILE *fp,*fps;
char *v,*v0;
char *mess;
{
   int i=2,n=1;
   int j;
   double d;

   if (readvar(varname,&j,&d,v,i,fp) >= 0) {
      if (mess[0]!='\0') fprintf(fps,"%s %s %s (readin)\n",mess,varname,v);
   }
   else {
      strcpy(v,v0); n=0; 
      if (mess[0]!='\0') fprintf(fps,"%s %s %s (default)\n",mess,varname,v);
   }
   return(n);
}
/*--------------------------------------------------------------------------*/
char *asplit(str,indstr,ind)
/* 
   Split string str in two parts on both sides of string "[indstr]";
   then place the numerical value ind in place of indstr. 
   Return the new string. 
   The routine looks like it does because I once thought the string functions
   caused tilts, so I removed almost all of them.
*/
char *str,*indstr;
int ind; 
{
   char index[10];
   char nindex[10];
   char b1[120],b2[120];
   static char buf[120];   /* This static is extremely important ! */
   int i,j,k;
   size_t len,len1,len2;
   char c;

   sprintf(index,"[%s]",indstr);
   sprintf(nindex,"[%d]",ind);
   len=strlen(str);
   len1=strlen(index);
   len2=strlen(nindex);

   for (i=0;i<len;i++) {
      if (sequal(str+i,index,len1)) break;
      b1[i]=str[i];
   }
   b1[i]='\0';
   if (i>=len) {
      fprintf(stderr,"asplit Error: index string not found: str %s index %s\n",str,index);
   }
   /*  printf("asplit1 %s %s %d %d %d %d\n",b1,str,i,len,len1,len2); oflush(); */
   for (j=0;j<len-i;j++) {
      b2[j]=str[i+len1+j];
      if (str[i+len1+j]=='\0' || str[i+len1+j]=='\n') break;
   }
   b2[j]='\0';
   for (k=0; k<=i+len2+j; k++) {
      if (k<i) c=str[k];
      if (k>=i && k<i+len2) c=nindex[k-i];
      if (k>=i+len2) c=b2[k-i-len2];
      buf[k]=c;
   }
   return(buf);

}
/*--------------------------------------------------------------------------*/

/*                                                                         */
/* This subroutine returns either an integer or a float in arg2 or arg3.   */
/* The file position will remain the old if any error is                   */
/* encountered, the new if the read is successful. The routine is case-    */
/* sensitive.                                                              */
/* The search always starts from the current file pointer position.        */
/* The search is through the entire file, starting                         */
/* anew from the beginning. Each new item must be separated by a space, tab*/
/* or newline from the previous. The name of each item may be followed by  */
/* any number of characters before the '=', which must be followed by the  */
/* value within 20 characters                                              */
/*
   flag:       0 float
               1 integer
               2 string
*/
/*
  Return values:
	>=0	Read successfull, position after read
	-1 	Name not found
	-2	'=' not found
	-3	File read error, or file doesn't exist
	-4	Name of variable too long
	-5	Value not found
*/
/*
   This routine was originally written for moldy.f, thus the horrifying
   code. For mdh.c, the flag=2 option for string reading was added.
   Maximum readin string length is 79 chars.
   Strings are read to the next space, tab or newline after the string,
   or as everything between two "".

*/
int readvar(varname,in,fl,str,flag,fp)
char *varname;
int *in;
double *fl;
char *str;
int flag;
FILE *fp;
{
   char buf[80],c,prevchar = 'a';
   int i,j,l=0,n,nsearched=2000000000;
   int originalposition,pos;
   
   static int firsttime=1;

   if (fp==NULL) return(-3);
   if (firsttime) {
      firsttime=0;
      fseek(fp,0,0);
   }

   *in=0; *fl=0.0;   /* Necessary to prevent tilts in some architectures */

   l=0; while (varname[l]!=':' && l<=30) l++;  l++;    /* Get length of varname */
   if (l>30) return(-4);          /* name too long */
   originalposition=ftell(fp);

/* Search for varname */

   pos=ftell(fp); 
   searchname(fp,varname,l,&n,&nsearched,&pos);
   if (n == -1) { fseek(fp,originalposition,0); return(-3); } /*file input error*/
   if (n<l && nsearched > 1000000000) {
      pos=0;    /* Search from beginning of file */
      searchname(fp,varname,l,&n,&nsearched,&pos);
   }
   if (n == -1) { fseek(fp,originalposition,0); return(-3); } /*file input error*/

   if (n==0) { fseek(fp,originalposition,0); return(-2); } /* name not found */

/* Name found, skip '=' and spaces */

   while ( (n=fread(&c,1,1,fp))>0 && (c != '=') ) fseek(fp,++pos,0);
   
   if (n == -1) { fseek(fp,originalposition,0); return(-3); } /*file input error*/
   if (n == 0)  { fseek(fp,originalposition,0); return(-1); } /* '=' not found  */
   if (!skipspaces(fp) && flag != 2) { 
      fseek(fp,originalposition,0); return(-5); /*value not found */
   }

/* Read in value */
   if (n<0) puts("readvar.c: This shouldn't be possible !"); 
  
   n=fread(buf,1,79,fp);   
   if (n == -1) {fseek(fp,originalposition,0); return(-3);} /*file input error*/
   if (n == 0)  {fseek(fp,originalposition,0); return(-5);} /*value not found */

   if (flag==1) sscanf (buf,"%d",in);            /* integer */
   else if (flag==2) {                           /* string */
      if (buf[0]!='"') sscanf(buf,"%s",str);
      else for(i=1;i<=79 && buf[i]!='"';i++) str[i-1]=buf[i];
   }
   else sscanf (buf,"%lf",fl);                   /* double */
   fseek(fp,-n+1,1);
   while( (n=fread(&c,1,1,fp))>0 && c!=' ' && c!='\n' && c!='\t');
   
   return(ftell(fp));
}
skipspaces(fp)
FILE *fp;
{
   char c;
   
   do { fread(&c,1,1,fp); }
   while (c == ' ' || c == '\n' || c == '\t');
   fseek(fp,-1,1);
   
   if ( (c>='0' && c<='9') || c=='.' || c=='+' || c=='-' ) return(1);
   return(0);
}

searchname(fp,varname,l,n,nsearched,pos)
/*
  Subroutine searches for name starting from position pos, until it
encounters the name, end of file or an error. Returns the position where the
name was found, or n<0 if it wasn't.
*/
FILE *fp;
char *varname;
int l,*pos,*n,*nsearched;
{
   unsigned char prevchar;
   int matches;
   unsigned char *c,buf[80];

   prevchar=' '; *n=1; matches=0;
   fseek(fp,*pos,0);
   while ((*nsearched)-- && matches<l) {
      matches=0;
      c=buf;
      while (matches<l) {
         if (*n=fread(c,1,1,fp) <= 0) return;       
         if (isspace(*c) || !isspace(prevchar)) break;
         if (*c==varname[matches]) { c++; matches++; }
         else break;
      }
      prevchar=buf[0];
      (*pos)++;
      fseek(fp,*pos,0);
   }
   if (matches==l) {
      *n=matches;
      *pos=*pos-1+l;
   }
   fseek(fp,*pos,0);
}
sequal(s1,s2,l)
/*
   Return 1 if strings s1 and s2 are equal to l places,
   else return 0.
*/
char s1[],s2[];
int l;
{
   while( (l-- >0) && *s1++ == *s2++) ;
   if (l== -1) return(1);
   return(0); 
}
