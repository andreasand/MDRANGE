
/* 
   Definitions for float (not double) math functions.
  
  These are retained only to enable calling the faster float math
  functions on some systems: in the default "double" version of
  the code they shouldn't be used of course !

*/

#define sqrtf(a)  sqrt((double) (a))
#define expf(a)   exp((double) (a))
#define sinf(a)   sin((double) (a))
#define cosf(a)   cos((double) (a))
#define tanf(a)   tan((double) (a))
#define asinf(a)  asin((double) (a))
#define acosf(a)  acos((double) (a))
#define atanf(a)  atan((double) (a))
#define logf(a)   log((double) (a))
#define powf(a,b) pow((double) (a),(double) (b))


