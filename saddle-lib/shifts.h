
typedef struct {
  unsigned int i;
  unsigned int j;
  double value;
} SparseRow;

void matrixVector(int n, SparseRow *a, int na, double *in, double *out);
void vectorMatrix(int n, SparseRow *a, int na, double *in, double *out);


/* FORTRAN to C data type matching */
typedef double doublereal;
typedef int integer;
typedef int fsub;
typedef fsub (*U_fp)(); /* Unknown procedure type (function name) */
typedef /* Subroutine */ fsub (*S_fp)();
typedef doublereal (*D_fp)();
typedef integer ftnlen;
typedef integer logical;
/* Used "nm minresqlpModule.o" to determine this */
#define MINRESQLP __minresqlpmodule_MOD_minresqlp

extern int MINRESQLP(
    integer *n, S_fp aprod, doublereal *b, doublereal *shift, S_fp msolve,
    logical *disable,  integer *nout,
    integer *itnlim, doublereal *rtol, doublereal *maxxnorm,
    doublereal *trancond, doublereal *acondLim,
    doublereal *x, integer *istop, integer *itn, doublereal *rnorm,
    doublereal *arnorm, doublereal *xnorm, doublereal *anorm,
    doublereal *acond);


void dynamicInit(unsigned int n, SparseRow *gauge,
		 unsigned int gaugeDimension, unsigned int gaugeElements,
		 integer itnlim, doublereal rtol);
void dynamicProject(unsigned int n, double *in, double *out);
int gaugeProduct(integer *vectorLength, doublereal *x, doublereal *y);

