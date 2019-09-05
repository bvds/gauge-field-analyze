#include <cjson/cJSON.h>

#ifdef USE_LIBRSB
#include <rsb.h>         /* for rsb_lib_init */
typedef struct rsb_mtx_t SparseRow;
#else
typedef struct {
    unsigned int i;
    unsigned int j;
    double value;
} SparseRow;
#endif

void matrixVector(const int n, const SparseRow *a, const int na,
                  const double *in, double *out);
void vectorMatrix(const int n, const SparseRow *a, const int na,
                  const double *in, double *out);


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
#define DGEMV dgemv_

extern int MINRESQLP(
    const integer *n, S_fp aprod, const doublereal *b,
    const doublereal *shift, S_fp msolve, S_fp userOrtho, const logical *disable,
    const integer *nout, const integer *itnlim, const doublereal *rtol,
    const doublereal *maxxnorm, const doublereal *trancond,
    const doublereal *acondLim,
    doublereal *x, integer *istop, integer *itn, doublereal *rnorm,
    doublereal *arnorm, doublereal *xnorm, doublereal *anorm,
    doublereal *acond);
// BLAS2 routine
extern int DGEMV(const char *, const integer *, const integer *,
                 const doublereal *, const doublereal *, const integer *,
                 const doublereal *, const integer *, const doublereal *,
                 doublereal *, const integer *);

void dynamicInit(unsigned int n, SparseRow *gauge,
		 unsigned int gaugeDimension, unsigned int gaugeElements,
		 cJSON *options);
void dynamicProject(const int n, double *in, double *out);
int gaugeProduct(const integer *vectorLength, const doublereal *x,
                 doublereal *y);
void printDynamicStats();


void hessInit(unsigned int n, SparseRow *hess, unsigned int hessElements);
void hessOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam);
void largeShifts(int n, double *initialVector, cJSON *options,
		 double **vals, double **vecs, unsigned int *nvals);
void testOp(const int n, double *grad);


void cutoffNullspace(unsigned int n, unsigned int nvals, cJSON *options,
                     double *grad,
                     double *vals, double *vecs, unsigned int *nLargeShifts);


void linearInit(unsigned int n, SparseRow *hess, int hessElements,
                double *vecs, int nvecs);
int hessProduct(integer *vectorLength, doublereal *x, doublereal *y);
void linearSolve(integer n, double *b, cJSON *options, double *x);
void userOrtho(char *action, integer *n, double *y);
void largeShiftProject(integer n, double *y);
