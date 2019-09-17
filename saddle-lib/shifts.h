#include <cjson/cJSON.h>

#ifdef USE_LIBRSB
#include <rsb.h>         /* for rsb_lib_init */
typedef struct rsb_mtx_t SparseMatrix;
#elif defined(USE_MKL)
#include "mkl.h"
typedef struct {
    sparse_matrix_t a;
    double *value;
    MKL_INT *row;
    MKL_INT *column;
    MKL_INT rows;
    MKL_INT columns;
    MKL_INT nonzeros;
} SparseMatrix;
#else
typedef struct {
    unsigned int i;
    unsigned int j;
    double value;
} SparseRow;

typedef struct {
    SparseRow *data;
    int nonzeros;
    int rows;
    int columns;
} SparseMatrix;
#endif

#ifdef USE_LIBRSB
int rows(SparseMatrix *matrix);
int columns(SparseMatrix *matrix);
int nonzeros(SparseMatrix *matrix);
#else
#define rows(A) A->rows
#define columns(A) A->columns
#define nonzeros(A) A->nonzeros
#endif


void matrixVector(const char type, const SparseMatrix *a,
                  const double *in, double *out);
void vectorMatrix(const char type, const SparseMatrix *a,
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
#define DNRM2 dnrm2_
#define DAXPY daxpy_

extern int MINRESQLP(
    const integer *n, S_fp aprod, const doublereal *b,
    const doublereal *shift, S_fp msolve, S_fp userOrtho, const logical *disable,
    const integer *nout, const integer *itnlim, const doublereal *rtol,
    const doublereal *abstol,
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
extern doublereal DNRM2(const integer *, const doublereal *, const integer *);
extern int DAXPY(const integer *n,  const doublereal *alpha,
                 const doublereal *x, const integer *incx,
                 const doublereal *y, const integer *incy); 

void dynamicInit(SparseMatrix *gauge, cJSON *options);
void dynamicProject(const int n, double *v, double *normDiff);
void gaugeOp(const int nrow, const int ncol, const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam);
int gaugeProduct(const integer *vectorLength, const doublereal *x,
                 doublereal *y);
void dynamicClose();


void hessOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam);
void hessOp2(const int nrow, const int ncol, const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam);
void largeShifts(SparseMatrix *, double *initialVector, cJSON *options,
		 double **vals, double **vecs, unsigned int *nvals);
void testOp(SparseMatrix *hess, double *grad);


void cutoffNullspace(unsigned int n, unsigned int nvals, cJSON *options,
                     double *grad,
                     double *vals, double *vecs, unsigned int *nLargeShifts);


void linearInit(SparseMatrix *hess, double *vecs, int nvecs);
int hessProduct(integer *vectorLength, doublereal *x, doublereal *y);
void linearSolve(integer n, double *b, cJSON *options, double *x);
void userOrtho(char *action, integer *n, double *y);
void largeShiftProject(integer n, double *y);
