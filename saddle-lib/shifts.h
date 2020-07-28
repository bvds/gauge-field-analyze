#include <cjson/cJSON.h>
#ifdef USE_MKL
#include "mkl.h"
#endif
#include "fortran.h"

#ifdef USE_MPI
#include <mpi.h>
#define _MPI_Comm MPI_Comm
#define TIME_TYPE double
#define SET_TIME(A) A=MPI_Wtime()
#define ADD_TIME(A,B,C) do{B=MPI_Wtime();A+=B-C;}while(0)
#else
#include <sys/time.h>
#define TIME_TYPE struct timeval
#define SET_TIME(A) gettimeofday(&A,NULL)
#define ADD_TIME(A,B,C) do{gettimeofday(&B,NULL);(A)+=tdiff(B,C);}while(0)
#define _MPI_Comm void *
#endif

/*
       Memory allocation
*/
#ifdef USE_MPI
/* MPI has its own versions of memory allocation,
   which may allow memory mapping for access. 

   Some of the vectors accessed by MPI in matrixVector() 
   come from external library routines:
   MINRESQLP(... hessProduct(... x, y) ... userOrtho(... y) ...)
   MINRESQLP(... gaugeProduct(... x ...) ...)
   trlan(hessOp2(... xin ... yout ...), dynamicProject(... v ...) ...)
*/
#define MALLOC(A, B) MPI_Alloc_mem(B,MPI_INFO_NULL,&(A))
#define FREE(A) MPI_Free_mem(A)
#else
#define MALLOC(A, B) A=malloc(B)
#define FREE(A) free(A)
#endif


/*
      Sparse matrix data structures
*/
// Index for sparse matrix rows/columns, nonzero elements
typedef unsigned int mat_int;
#ifdef USE_MPI
// For MPI calls, set matching type
#define _MPI_MAT_INT MPI_UNSIGNED
#endif
#ifdef USE_TASK
    typedef struct {
        mat_int start;
        mat_int end;
        int doRank;
        int sendTo;
        int receiveFrom;
        mat_int receiveSize;
    } TaskList;
#endif

/*  The more intuitive construct is actually undefined;
    see http://lists.llvm.org/pipermail/cfe-commits/Week-of-Mon-20160118/147239.html
*/
#if defined(USE_MKL) && !defined(USE_MPI) && !defined(USE_BLOCK) && !defined(USE_TASK)
#define USE_MKL_MATRIX 1
#else
#define USE_MKL_MATRIX 0
#endif

typedef struct {
    mat_int blocks;
    mat_int blockSize;
    mat_int rows;       // global number of rows
    mat_int columns;    // global number of columns
    mat_int rowParts;   // row partitions
    mat_int colParts;   // column partitions
    double *value;
#if USE_MKL_MATRIX
    MKL_INT *i;
    MKL_INT *j;
    struct matrix_descr descr;
    sparse_matrix_t a;
#else
    mat_int *i;
    mat_int *j;
    mat_int *ip;        // beginning and end of a row
    mat_int iCount;     // local number of distinct rows
#endif
#ifdef USE_MPI
    MPI_Comm mpicom;
#endif
    double mpiTime;
    double localTime;
    int count;
#ifdef USE_TASK
    TaskList* task;
    int taskCount;
    double *gather[2];
#else
    double *gather;
    mat_int lowerColumn;
#endif
} SparseMatrix;


/*
            shifts.c
*/
int getPrintDetails(cJSON *jopts, int defaultLevel);
#ifndef USE_MPI
double tdiff(struct timeval t1, struct timeval t0);
#endif


/*
            dynamic.c
 */
void dynamicInit(const mat_int nrow, const mat_int ncol,
                 SparseMatrix *gauge, SparseMatrix *gaugeT,
                 cJSON *options, _MPI_Comm mpicom);
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
		 double **eval, double **evec, int *nvals, _MPI_Comm mpicom);
void hessOp(const int nrow, const int ncol,
            const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam);
void hessOp2(const int nrow, const int ncol,
             const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam);
void testOp(SparseMatrix *hess, const mat_int nrow,
            double *grad, _MPI_Comm mpicom);


/*
         cutoff.c
*/
void cutoffNullspace(mat_int n, int nvals, cJSON *options,
                     double *grad,
                     double **vals, double **vecs, int *nLargeShifts,
                     cJSON *jout, _MPI_Comm mpicom);


/*
            linear.c
 */
void linearInit(SparseMatrix *hess, const mat_int nrow, double *vecs,
                int nvecs, _MPI_Comm mpicom);
void hessProduct(integer *n, doublereal *x, doublereal *y);
void linearSolve(const integer n, double *b, cJSON *options, double *x);
void userOrtho(char *action, integer *n, double *y);
void largeShiftProject(integer n, double *y);


/*
         sort.c
*/
void sortMatrixChunks(SparseMatrix *mat, const mat_int chunk);
void sortMatrixLocal(SparseMatrix *mat, int wrank, int wsize, int debug);


/*
         matrix.c
*/
int indexRank(const mat_int i, const int wsize,
              const mat_int n, const mat_int partitions);
mat_int localSize(const unsigned int wrank, const int wsize,
                  const mat_int n, const mat_int partitions);
mat_int rankIndex(const unsigned int wrank, const int wsize,
                  const mat_int n, const mat_int partitions);
mat_int maxLocalSize(const int wsize,
                     const mat_int n, const mat_int partitions);
void rankSanityTest(const mat_int n, const mat_int partitions);
void sparseMatrixRead(SparseMatrix *mat, char *fileName, char descr,
                      int tFlag, int blockSize, mat_int rowParts,
                      mat_int colParts,
                      int chunkSize, int debug, _MPI_Comm mpicom);
void sparseMatrixStats(SparseMatrix *mat, char *matrixName);
void sparseMatrixFree(SparseMatrix *mat);
void matrixVector(SparseMatrix *a,
                  const mat_int lin, const doublereal *in,
                  const mat_int lout, doublereal *out);
void testMatrixVector(SparseMatrix *mat, double *in);
