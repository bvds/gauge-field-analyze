#include <cjson/cJSON.h>

// Index sparse matrix rows/columns
typedef unsigned int mat_int;

#ifdef USE_MKL
#include "mkl.h"
#define MALLOC_ALIGN 64
#endif
#ifdef USE_BLOCK
typedef struct {
    double *value;
    mat_int *i;
    mat_int *j;
    mat_int blocks;
    mat_int blockSize;
    mat_int nonzeros;
    mat_int rows;
    mat_int columns;
    char descr;
} SparseMatrix;
#elif defined(USE_MKL)
typedef struct {
    sparse_matrix_t a;
    struct matrix_descr descr;
    double *value;
    mat_int blockSize;
    MKL_INT *i;
    MKL_INT *j;
    mat_int rows;
    mat_int columns;
    mat_int nonzeros;
} SparseMatrix;
#else
typedef struct {
    double *value;
    mat_int *i;
    mat_int *j;
    char descr;
    mat_int nonzeros;
    mat_int rows;
    mat_int columns;
} SparseMatrix;
#endif

#define rows(A) A->rows
#define columns(A) A->columns
#define nonzeros(A) A->nonzeros


/* FORTRAN to C data type matching */
typedef double doublereal;
typedef int integer;
typedef void fsub;  // gfortran convension
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
void dynamicProject(const integer n, double *v, double *normDiff);
void gaugeOp(const int nrow, const int ncol, const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam);
void gaugeProduct(const integer *vectorLength, const doublereal *x,
                 doublereal *y);
void dynamicClose();
void matrixVector(const SparseMatrix *a,
                  const doublereal *in, doublereal *out);
void vectorMatrix(const SparseMatrix *a,
                  const doublereal *in, doublereal *out);


void hessOp(const int nrow, const int ncol,
            const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam);
void hessOp2(const int nrow, const int ncol,
             const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam);
void eigenInit(SparseMatrix *hess);
void largeShiftsCheckpoint(double *initialVector, cJSON *options,
		 double **vals, double **vecs, int *nvals);
void largeShifts(double *initialVector, cJSON *options,
		 double **vals, double **vecs, int *nvals);
void testOp(SparseMatrix *hess, double *grad);


void cutoffNullspace(mat_int n, mat_int nvals, cJSON *options,
                     double *grad,
                     double **vals, double **vecs, mat_int *nLargeShifts);


void linearInit(SparseMatrix *hess, double *vecs, int nvecs);
void hessProduct(integer *vectorLength, doublereal *x, doublereal *y);
void linearSolve(integer n, double *b, cJSON *options, double *x);
void userOrtho(char *action, integer *n, double *y);
void largeShiftProject(integer n, double *y);

void sortMatrix(SparseMatrix *mat, const mat_int chunk);
