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

void eigenInit(SparseMatrix *hess) {
    eigenData.matrix = hess;
}

/*
   Wrapper to optionally read/write checkpoint file.
   For debugging.
*/

void largeShiftsCheckpoint(double *initialVector, cJSON *options,
                           double **vals, double **vecs, int *nvals) {
    cJSON *tmp;
    char *checkpoint;
    FILE *fp;
    int k;
    mat_int i, nn, n = rows(eigenData.matrix);
    tmp = cJSON_GetObjectItemCaseSensitive(options, "readCheckpoint");
    if(tmp == NULL) {
        largeShifts(initialVector, options, vals, vecs, nvals);
        tmp = cJSON_GetObjectItemCaseSensitive(options, "writeCheckpoint");
        if(tmp != NULL) {
            checkpoint = tmp->valuestring;
            printf("Opening output file %s\n", checkpoint);
            if((fp = fopen(checkpoint, "w")) == NULL) {
                fprintf(stderr, "Can't open file %s\n", checkpoint);
                exit(555);
            }
            fprintf(fp,"%u %u\n", n, *nvals);
            for(i=0; i<n; i++){
                fprintf(fp, "%.17e\n", (*vals)[i]);
            }
            for(i=0; i<*nvals*n; i++){
                fprintf(fp, "%.17e\n", (*vecs)[i]);
            }
            fclose(fp);
        }
    } else {
        checkpoint = tmp->valuestring;
        printf("Opening input file %s\n", checkpoint);
        if((fp = fopen(checkpoint, "r")) == NULL) {
            fprintf(stderr, "Can't open file %s\n", checkpoint);
            exit(555);
        }
        k = fscanf(fp, "%u%d", &nn, nvals);
        assert(nn==n);
        assert(k==2);
        *vals = malloc(n*sizeof(double));
        *vecs = malloc(n*(*nvals)*sizeof(double));
        for(i=0; i<n; i++) {
            k = fscanf(fp, "%le", (*vals)+i);
            assert(k==1);
        }
        for(i=0; i<n*(*nvals); i++) {
            k = fscanf(fp, "%le", (*vecs)+i);
            assert(k==1);
        }
        fclose(fp);
    }
}

/*
  Returns the absolute values of the lowest eigenvalues of
  the Hessian.
 */
void largeShifts(double *initialVector, cJSON *options,
		 double **eval, double **evec, int *nvals) {
#ifdef USE_PRIMME
    primme_params primme;
    double *rnorm;   /* Array with the computed eigenpairs residual norms */
    double targetShifts[1];
#else
    // Let trlan figure out the work array allocation.
    const int lwrk = 0;
    const int iguess = 1; // Use supplies the intial vector in evec
    double *wrk = NULL;
    int mev, maxlan, maxmv;
    double tol = -1.0; // Use default tolerance: sqrt(machine epsilon)
    /* Strategies 3 & 4 are for matrices where the matrix-vector 
       multiply is relatively expensive.  */
    int restart = 4;
    trl_info info;
    int k;
#endif
    SparseMatrix *hess = eigenData.matrix;
    int nrow = rows(hess); // number of rows on this processor
    int ned = 1;  // number of requested eigenpairs.
    int lohi = -1; // lowest or highest abs value eigenpairs.
    int i, ierr, printDetails = 1;
    cJSON *tmp;
    clock_t t1;
    time_t t2, tf;

    t1 = clock();
    time(&t2);
    
    eigenData.ops = 0;
    eigenData.sumNcol = 0;
    eigenData.maxNcol = 0;
    eigenData.tcpu_mv = 0;

    assert(columns(hess) == rows(hess));
    assert(eval != NULL);
    assert(evec != NULL);

    // Would need a separate flag for the lohi == 0 case.
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenPairs");
    if(cJSON_IsNumber(tmp)) {
        ned = abs(tmp->valueint);
        lohi = tmp->valueint<0?-1:1;
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "printDetails");
    if(cJSON_IsNumber(tmp)) {
        printDetails = tmp->valueint;
    } else if(cJSON_IsBool(tmp)) {
        printDetails = cJSON_IsTrue(tmp)?1:0;
    }    
#ifdef USE_PRIMME
    /* Set default values in PRIMME configuration struct */
    primme_initialize(&primme);
    primme.matrixMatvec = hessOpPrimme;
    primme.n = rows(hess);
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
        maxmv = rows(hess) * ned;
    }
    /* maxLanczosVecs is maximum number of Lanczos vectors to use.
       This sets an upper bound on the memory used. */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxLanczosVecs");
    if(cJSON_IsNumber(tmp)) {
        maxlan = tmp->valueint;
    } else {
        maxlan = rows(hess);  // Since we do reorthogonalization.
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "restartStrategy");
    if(cJSON_IsNumber(tmp)) {
        restart = tmp->valuedouble;
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        tol = tmp->valuedouble;
    }

    // Used by HessOp2
#ifdef USE_MKL
    eigenData.z = mkl_malloc(rows(hess) * sizeof(double), MALLOC_ALIGN);
#else
    eigenData.z = malloc(rows(hess) * sizeof(double));
#endif
    mev = ned; // Allocate memory for the number of requested eigenpairs
    *eval = (double *) malloc(mev*sizeof(double));
    *evec = (double *) malloc(mev*nrow*sizeof(double));
    trl_init_info(&info, nrow, maxlan, lohi, ned, tol, restart, maxmv, 0);
    trl_set_iguess(&info, 0, iguess, 0, NULL);
    memset(*eval, 0, mev*sizeof(double));
    // Provide the initial vector
    for( i=0; i<nrow; i++ ) {
        (*evec)[i] = initialVector==NULL?1.0:initialVector[i];
    }

    // call TRLAN to compute the eigenvalues
    trlan(hessOp2, dynamicProject, &info, nrow, mev, *eval, *evec,
          nrow, lwrk, wrk);
    *nvals = info.nec;
    ierr = info.stat;
    if(printDetails > 1) {
        /* This estimate of flops doesn't include dynamicProject() call. 
           Assuming 1 addition and 1 multiply per nonzero element,
           2 matrix multiplies. */
        trl_print_info(&info, 4*nonzeros(hess));
    } else if(printDetails > 0) {
        trl_terse_info(&info, stdout);
    }

#ifdef USE_MKL
    mkl_free(eigenData.z);
#else
    free(eigenData.z);
#endif
#endif

    time(&tf);
    if(printDetails > 0) {
        printf("largeShifts: ops=%i, sumNcol=%i, maxNcol=%i\n",
               eigenData.ops, eigenData.sumNcol, eigenData.maxNcol);
        printf("largeShifts: %i matrix-vector ops, "
               "%i reorthogonalizations in %.2f sec (%li wall)\n",
#ifdef USE_PRIMME
               0, 0,
#else
               info.matvec, info.north,
#endif
               (clock()-t1)/(float) CLOCKS_PER_SEC, tf-t2);
        printf("largeShifts:  %.2fmv cpu sec for ops\n",
               eigenData.tcpu_mv/(float) CLOCKS_PER_SEC);
        fflush(stdout);
    }

    if(ierr < 0 || *nvals < ned) {
        fprintf(stderr, "large shift exit with stat=%i, "
                "finding %i of %i eigenpairs.\n",
               ierr, *nvals, ned);
        exit(8);
    }

#ifndef USE_PRIMME
    /* We have calculated eigenvalues for hess^2.
       Convert to abs(eigenvalues) of hess itself. */
    for(k=0; k<*nvals; k++) {
        (*eval)[k] = sqrt((*eval)[k]);
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
    assert(columns(eigenData.matrix) == abs(nrow));
    assert(rows(eigenData.matrix) == abs(nrow));
    assert(mvparam == NULL);
    int k;
    clock_t t0;

    eigenData.ops += 1;
    eigenData.sumNcol += ncol;
    if(eigenData.maxNcol<ncol)
        eigenData.maxNcol = ncol;

    for(k=0; k<ncol; k++) {
        t0 = clock();
        matrixVector(eigenData.matrix, xin+k*ldx, yout+k*ldy);
        eigenData.tcpu_mv += clock() - t0;
        // Use the fact that "in" and "out" can overlap.
        dynamicProject(nrow, yout+k*ldy, NULL);
    }
}

/* The extra parameter mvparam is not used in this case. */
void hessOp2(const int nrow, const int ncol,
             const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam) {
    assert(columns(eigenData.matrix) == abs(nrow));
    assert(rows(eigenData.matrix) == abs(nrow));
    assert(mvparam == NULL);
    int k;
    clock_t t0;

    eigenData.ops += 1;
    eigenData.sumNcol += ncol;
    if(eigenData.maxNcol<ncol)
        eigenData.maxNcol = ncol;

    for(k=0; k<ncol; k++) {
        t0 = clock();
        matrixVector(eigenData.matrix, xin+k*ldx, eigenData.z);
        eigenData.tcpu_mv += clock() - t0;
        // Use the fact that "in" and "out" can overlap.
        dynamicProject(nrow, eigenData.z, NULL);
        // Apply a second time
        t0 = clock();
        matrixVector(eigenData.matrix, eigenData.z, yout+k*ldy);
        eigenData.tcpu_mv += clock() - t0;
        // Use the fact that "in" and "out" can overlap.
        dynamicProject(nrow, yout+k*ldy, NULL);
    }
}

// Debug print of dynamic part of hess.grad
void testOp(SparseMatrix *hess, double *grad) {
    double *y;
    mat_int i;

    y = malloc(rows(hess) * sizeof(double));
    hessOp(rows(hess), 1, grad, 1, y, 1, NULL);
    printf("Dynamic part of hess.grad\n");
    for(i=0; i<rows(hess); i++) {
        printf("  %le\n", y[i]);
    }
    free(y);
}
