#include <cjson/cJSON.h>
#ifdef USE_MKL
#include "mkl.h"
#endif
#include "fortran.h"

/*
       Memory allocation
*/
#ifdef USE_MPI
#include <mpi.h>
/* TODO:  MPI has its own versions of memory allocation,
   which may be more efficient for shared memory windows. */
#define MALLOC(A) malloc(A)
#define REALLOC(A, B) realloc(A, B)
#define FREE(A) free(A)
#elif defined(USE_MKL)
#define MALLOC_ALIGN 64
#define MALLOC(A) mkl_malloc(A, MALLOC_ALIGN)
#define REALLOC(A, B) mkl_realloc(A, B)
#define FREE(A) mkl_free(A)
#else
#define MALLOC(A) malloc(A)
#define REALLOC(A, B) realloc(A, B)
#define FREE(A) free(A)
#endif


/*
      Sparse matrix data structures
*/
// Index for sparse matrix rows/columns, nonzeros
typedef unsigned int mat_int;
#ifdef USE_MPI
// For MPI calls, set matching type
#define _MPI_MAT_INT MPI_UNSIGNED_INT
#endif

typedef struct {
    mat_int nonzeros;
    mat_int rows;
    mat_int columns;
    double *value;
#ifdef USE_BLOCK
    mat_int *i;
    mat_int *j;
    char descr;
    mat_int blocks;
    mat_int blockSize;
#elif defined(USE_MKL) && !defined(USE_MPI)
    MKL_INT *i;
    MKL_INT *j;
    struct matrix_descr descr;
    sparse_matrix_t a;
    mat_int blockSize;
#else
    mat_int *i;
    mat_int *j;
    char descr;
#endif
} SparseMatrix;


/*
            dynamic.c
 */
void dynamicInit(const mat_int nrow, const mat_int ncol,
                 SparseMatrix *gauge, cJSON *options, void *mpicomp);
void dynamicProject(const integer n, double *v, double *normDiff);
void gaugeOp(const int nrow, const int ncol, const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam);
void gaugeProduct(const integer *vectorLength, const doublereal *x,
                 doublereal *y);
void dynamicClose();


/*
              eigensystem.c
 */
void largeShifts(SparseMatrix *hess, cJSON *options,
                 const mat_int nrow, double *initialVector,
		 double **eval, double **evec, int *nvals, void *mpicomp);
void hessOp(const int nrow, const int ncol,
            const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam);
void hessOp2(const int nrow, const int ncol,
             const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam);
void testOp(SparseMatrix *hess, const mat_int nrow,
            double *grad, void *mpicomp);


/*
         cutoff.c
*/
void cutoffNullspace(mat_int n, int nvals, cJSON *options,
                     double *grad,
                     double **vals, double **vecs, int *nLargeShifts,
                     void *mpicom);


/*
            linear.c
 */
void linearInit(SparseMatrix *hess, const mat_int nrow, double *vecs,
                int nvecs, void *mpicomp);
void hessProduct(integer *n, doublereal *x, doublereal *y);
void linearSolve(const integer n, double *b, cJSON *options, double *x);
void userOrtho(char *action, integer *n, double *y);
void largeShiftProject(integer n, double *y);


/*
         sort.c
*/
void sortMatrix(SparseMatrix *mat, const mat_int chunk);


/*
         matrix.c
*/
int indexRank(const mat_int i, const int wsize, const mat_int n);
mat_int localSize(const unsigned int wrank, const int wsize, const mat_int n);
mat_int rankIndex(const unsigned int wrank, const int wsize, const mat_int n);
void rankSanityTest(mat_int n);
void sparseMatrixRead(FILE *fp, SparseMatrix *mat, char *fileName, int wrank);
void sparseMatrixFree(SparseMatrix *mat);
void matrixVector(const SparseMatrix *a,
                  const mat_int lin, const doublereal *in,
                  const mat_int lout, doublereal *out, void *mpicomp);
void vectorMatrix(const SparseMatrix *a,
                  const mat_int lin, const doublereal *in,
                  const mat_int lout, doublereal *out, void *mpicomp);
