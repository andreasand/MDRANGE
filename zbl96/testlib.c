
#include "zbl96lib.h"
#include "testlib.h"

int main()
{
   double *S,x;
   unsigned int flag;
   int i;
   
   flag = DEFAULT;
      
   S = zbl96(14,14,0,28.086,2.3212,0,10,1,flag);

   for(x=0,i=0;x<=10;x+=1,i++)
      printf("%12.5f %12.5f\n",x,S[i]);

   free(S);
      
   exit(0);
}
