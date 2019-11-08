/*
  Find eigenvalues associated with large shifts.

  These would be among the eigenvalues with the smallest
  magnitude.  Thus, we calculate the smallest eigenvalues
  of the Hessian squared using thick-restart Lanczos.
  This is the largest time user in the calculation.

  One could use the Jacobi-Davidson algorithm
  on the Hessian itself.
  For instance, 
    PRIMME [Stathopoulos, 2007] C library
    https://github.com/primme/primme or

    JADAMILU [Bollhöfer and Notay, 2007]
    Does not have source code available.

  or use filtered methods

    The FEAST eigensolver.
    http://www.ecs.umass.edu/~polizzi/feast/

    EVSL
    https://github.com/eigs/EVSL

  Maybe going to block lanczos would help?  The fact that 
  the Hessian itself has nc^2-1 size blocks is not relevant.

    KSHELL
    https://sites.google.com/a/cns.s.u-tokyo.ac.jp/kshell/

  Finally, one could use "shift and invert" with the 
  shift = 0 to calculate the eigensystem of A^-1.
  However, we already know that solving the final linear
  system takes many Lanczos iterations.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "shifts.h"
#include "trlan.h"
#include "trl_comm_i.h"
#ifdef USE_PRIMME
#include "primme.h"
#endif

#ifdef USE_PRIMME
void hessOpPrimme(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
                  int *blockSize, primme_params *primme, int *ierr);
#endif
struct {
    SparseMatrix *matrix;
    double *z;
    int ops;
    int sumNcol;
    int maxNcol;
    clock_t tcpu_mv;
} eigenData;


/*
  Returns the absolute values of the lowest eigenvalues of
  the Hessian.
 */
void largeShifts(SparseMatrix *hess, cJSON *options,
                 const mat_int nrow, double *initialVector,
		 double **eval, double **evec, int *nvals,
                 _MPI_Comm mpicom) {
#ifdef USE_PRIMME
    primme_params primme;
    double *rnorm;   /* Array with the computed eigenpairs residual norms */
    double targetShifts[1];
#else
    // Let trlan figure out the work array allocation.
    const int lwrk = 0;
    const int iguess = 1; // Use supplies the intial vector in evec
    double *wrk = NULL;
    int mev, maxlan, maxmv, wrank, restart;
    double tol;
    trl_info info;
    mat_int k;
#endif
    int ned = 1;  // number of requested eigenpairs.
    int lohi = -1; // lowest or highest abs value eigenpairs.
    int i, ierr, printDetails;
    mat_int nonzeros = hess->blocks*hess->blockSize*hess->blockSize;
#ifdef USE_MPI
    MPI_Comm *mpicomp = &mpicom;
    MPI_Comm_rank(mpicom, &wrank);
#else
    assert(mpicom == NULL);
    void *mpicomp = NULL;
    wrank = 0;
#endif
    cJSON *tmp;
    clock_t t1, t1f;
    time_t t2, tf;

    t1 = clock();
    time(&t2);
    
    eigenData.matrix = hess;
    eigenData.ops = 0;
    eigenData.sumNcol = 0;
    eigenData.maxNcol = 0;
    eigenData.tcpu_mv = 0;

    assert(eval != NULL);
    assert(evec != NULL);

    // Would need a separate flag for the lohi == 0 case.
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenPairs");
    if(cJSON_IsNumber(tmp)) {
        ned = abs(tmp->valueint);
        lohi = tmp->valueint<0?-1:1;
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "printDetails");
    if(cJSON_IsNumber(tmp))
        printDetails = tmp->valueint;
    else if(cJSON_IsBool(tmp))
        printDetails = cJSON_IsTrue(tmp)?1:0;
    else
        printDetails = 0;
#ifdef USE_PRIMME
    /* Set default values in PRIMME configuration struct */
    primme_initialize(&primme);
    primme.matrixMatvec = hessOpPrimme;
    primme.n = nrow;
    primme.numEvals = ned;
    primme.target = lohi<0?primme_closest_abs:primme_largest_abs;
    primme.numTargetShifts = 1;
    targetShifts[0] = 0.0;  // Eigenvalues near zero
    primme.targetShifts = targetShifts;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "rTolerance");
    primme.eps = cJSON_IsNumber(tmp)?abs(tmp->valuedouble):DBL_EPSILON;
    tmp = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    if (cJSON_IsNumber(tmp))
        primme.maxMatvecs = tmp->valueint;
    /* Set method to solve the problem */
    primme_set_method(PRIMME_DYNAMIC, &primme);
    /* DYNAMIC uses a runtime heuristic to choose the fastest method between
       PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
       set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */
    primme.initSize = 1; // Number of initial vectors to start with.
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxLanczosVecs");
    if(cJSON_IsNumber(tmp))
        primme.maxBasisSize = tmp->valueint;

    *eval = (double *) malloc(primme.numEvals*sizeof(double));
    *evec = (double *) malloc(primme.numEvals*primme.n*sizeof(double));
    rnorm = (double *) malloc(primme.numEvals*sizeof(double));
    // Provide the initial vector
    for( i=0; i<nrow; i++ ) {
        (*evec)[i] = initialVector==NULL?1.0:initialVector[i];
    }
    
    ierr = dprimme(*eval, *evec, rnorm, &primme);
    *nvals = primme.initSize;
    
    free(rnorm);
    primme_free(&primme);
#else

    /* We interpret maxIterations to be number of matrix-vector multiples.
       This sets an upper limit on time used. */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    if(cJSON_IsNumber(tmp)) {
        maxmv = tmp->valueint;
    } else {
        // Documented default (maxmv<0) seems to be broken?
        maxmv = hess->rows * ned;
    }
    /* maxLanczosVecs is maximum number of Lanczos vectors to use.
       This sets an upper bound on the memory used. */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxLanczosVecs");
    if(cJSON_IsNumber(tmp)) {
        maxlan = tmp->valueint;
    } else {
        maxlan = hess->rows;  // Since we do reorthogonalization.
    }
    /* Strategies 3 & 4 are for matrices where the matrix-vector 
       multiply is relatively expensive.  
       Benchmark test for 16^3, nc=3 lattice show 7 is fastest, with
       roughly 50% going to restarts and reorthogonalization
       and 50% to matrixVector().  */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "restartStrategy");
    restart = cJSON_IsNumber(tmp)?tmp->valuedouble:7;
    // Default is to use the default tolerance: sqrt(machine epsilon)
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "rTolerance");
    tol = cJSON_IsNumber(tmp)?tmp->valuedouble:-1.0; 

    // Used by HessOp2
    eigenData.z = MALLOC(nrow * sizeof(double));
    mev = ned; // Allocate memory for the number of requested eigenpairs
    *eval = (double *) malloc(mev*sizeof(double));
    *evec = (double *) malloc(mev*nrow*sizeof(double));
    trl_init_info(&info, nrow, maxlan, lohi, ned, tol, restart, maxmv,
                  mpicomp);
    if(wrank != 0)
        trl_set_debug(&info, 0, "eigensytem-trlan-");
    trl_set_iguess(&info, 0, iguess, 0, NULL);
    memset(*eval, 0, mev*sizeof(double));
    // Provide the initial vector
    for (k=0; k<nrow; k++ ) {
        (*evec)[k] = initialVector==NULL?1.0:initialVector[k];
    }

    // call TRLAN to compute the eigenvalues
    trlan(hessOp2, dynamicProject, &info, nrow, mev, *eval, *evec,
          nrow, lwrk, wrk);
    *nvals = info.nec;
    ierr = info.stat;
    if(printDetails > 1) {
        /* This estimate of flops doesn't include the dynamicProject() call.
           2 matrix multiplies. + and * as FLOPS. */
        trl_print_info(&info, 4*nonzeros);
    } else if(printDetails > 0 && wrank==0) {
        trl_terse_info(&info, stdout);
    }

    FREE(eigenData.z);
#endif

    time(&tf);
    t1f = clock() - t1;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &eigenData.tcpu_mv, 1, MPI_LONG,
                  MPI_SUM, mpicom);
    MPI_Allreduce(MPI_IN_PLACE, &t1f, 1, MPI_LONG,
                  MPI_SUM, mpicom);
#endif
    if(printDetails > 0 && wrank==0) {
        printf("largeShifts: ops=%i, sumNcol=%i, maxNcol=%i\n",
               eigenData.ops, eigenData.sumNcol, eigenData.maxNcol);
        printf("largeShifts: %i matrix-vector ops, "
               "%i reorthogonalizations in %.2f sec (%li wall)\n",
#ifdef USE_PRIMME
               0, 0,
#else
               info.matvec, info.north,
#endif
               t1f/(float) CLOCKS_PER_SEC, tf-t2);
        printf("largeShifts:  %.2fmv cpu sec for ops\n",
               eigenData.tcpu_mv/(float) CLOCKS_PER_SEC);
        fflush(stdout);
    }

    if(ierr < 0 || *nvals < ned) {
        fprintf(stderr, "%i:  large shift exit with stat=%i, "
                "finding %i of %i eigenpairs.\n",
                wrank, ierr, *nvals, ned);
        exit(8);
    }

#ifndef USE_PRIMME
    /* We have calculated eigenvalues for hess^2.
       Convert to abs(eigenvalues) of hess itself. */
    for(i=0; i<*nvals; i++) {
        (*eval)[i] = sqrt((*eval)[i]);
    }
#endif
}

#ifdef USE_PRIMME
void hessOpPrimme(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
                  int *blockSize, primme_params *primme, int *ierr) {
    (void)(blockSize);
    (void)(ierr);
    hessOp(primme->n, 1, x, *ldx, y, *ldy, NULL);
}
#endif

/* The extra parameter mvparam is not used in this case. */
void hessOp(const int nrow, const int ncol,
            const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(mvparam == NULL);
    int k;
    clock_t t0;

    eigenData.ops += 1;
    eigenData.sumNcol += ncol;
    if(eigenData.maxNcol<ncol)
        eigenData.maxNcol = ncol;

    for(k=0; k<ncol; k++) {
        t0 = clock();
        matrixVector(eigenData.matrix, nrow, xin+k*ldx, nrow, yout+k*ldy);
        eigenData.tcpu_mv += clock() - t0;
        dynamicProject(nrow, yout+k*ldy, NULL);
    }
}

/* The extra parameter mvparam is not used in this case. */
void hessOp2(const int nrow, const int ncol,
             const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam) {
    assert(mvparam == NULL);
    int k;
    clock_t t0;

    eigenData.ops += 1;
    eigenData.sumNcol += ncol;
    if(eigenData.maxNcol<ncol)
        eigenData.maxNcol = ncol;

    for(k=0; k<ncol; k++) {
        t0 = clock();
        matrixVector(eigenData.matrix, nrow, xin+k*ldx, nrow, eigenData.z);
        eigenData.tcpu_mv += clock() - t0;
        dynamicProject(nrow, eigenData.z, NULL);
        // Apply a second time
        t0 = clock();
        matrixVector(eigenData.matrix, nrow, eigenData.z, nrow, yout+k*ldy);
        eigenData.tcpu_mv += clock() - t0;
        // Use the fact that "in" and "out" can overlap.
        dynamicProject(nrow, yout+k*ldy, NULL);
    }
}

// Debug print of dynamic part of hess.grad
void testOp(SparseMatrix *hess, const mat_int nrow,
            double *grad, _MPI_Comm mpicom) {
    double *y;
    mat_int k;
    int wrank, wsize, i;
#ifdef USE_MPI
    MPI_Comm_size(mpicom, &wsize);
    MPI_Comm_rank(mpicom, &wrank);
#else
    assert(mpicom == NULL);
    wsize = 1;
    wrank = 0;
#endif

    eigenData.matrix = hess;

    y = malloc(nrow * sizeof(double));
    hessOp(nrow, 1, grad, 1, y, 1, NULL);
    printf("Dynamic part of hess.grad\n");
    for(i=0; i<wsize; i++) {
        if(i == wrank) {
            for(k = 0; k<nrow; k++) {
                printf("  %le (%i:  %u)\n", y[k], i, k);
            }
        }
#ifdef USE_MPI
        MPI_Barrier(mpicom);
#endif
    }
    free(y);
}
