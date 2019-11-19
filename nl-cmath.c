#include "newlisp.h"
#include "protos.h"

/* turn on for extra debugging output */
/* #define BAYES_DEBUG */
/* #define DEBUG */


#ifdef HAVE_COMPLEX
#include <complex.h>
typedef complex double complexreal;
#define Creal(z)  (creal(z))
#define Cimag(z)  (cimag(z))
#define Cnew(r,i) ((r)+Ci*(i))
#define Cadd(x,y) ((x)+(y))
#define Csub(x,y) ((x)-(y))
#define Cmul(x,y) ((x)*(y))
#define Cdiv(x,y) ((x)/(y))
#define Cconj(z)  (conj(z))
#endif

