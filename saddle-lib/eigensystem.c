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
                           double **vals, double **vecs, unsigned int *nvals) {
    cJSON *tmp;
    char *checkpoint;
    FILE *fp;
    int k;
    unsigned int i, nn, n = rows(eigenData.matrix);
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
        k = fscanf(fp, "%u%u", &nn, nvals);
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
		 double **eval, double **evec, unsigned int *nvals) {
    // Let trlan figure out the work array allocation.
    const int lwrk = 0;
    const int iguess = 1; // Use supplies the intial vector in evec
    double *wrk = NULL;
    int mev, maxlan, lohi, ned, maxmv;
    double tol = -1.0; // Use default tolerance: sqrt(machine epsilon)
    SparseMatrix *hess = eigenData.matrix;
    int nrow = rows(hess); // number of rows on this processor
    trl_info info;
    int i, printDetails = 1;
    /* Strategies 3 & 4 are for matrices where the matrix-vector 
       multiply is relatively expensive.  */
    int restart = 4;
    unsigned int k;
    cJSON *tmp;
    clock_t t1;
    time_t t2, tf;

    t1 = clock();
    time(&t2);

    // Used by HessOp
#ifdef USE_MKL
    eigenData.z = mkl_malloc(rows(hess) * sizeof(double), MALLOC_ALIGN);
#else
    eigenData.z = malloc(rows(hess) * sizeof(double));
#endif
    eigenData.ops = 0;
    eigenData.sumNcol = 0;
    eigenData.maxNcol = 0;
    eigenData.tcpu_mv = 0;

    // Would need a separate flag for the lohi == 0 case.
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenPairs");
    if(cJSON_IsNumber(tmp)) {
        ned = abs(tmp->valueint);
        lohi = tmp->valueint<0?-1:1;
    } else {
        ned = 1;
        lohi = -1;
    }
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
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "printDetails");
    if(cJSON_IsNumber(tmp)) {
        printDetails = tmp->valueint;
    } else if(cJSON_IsBool(tmp)) {
        printDetails = cJSON_IsTrue(tmp)?1:0;
    }

    assert(columns(hess) == rows(hess));
    assert(eval != NULL);
    assert(evec != NULL);

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

    time(&tf);
    if(printDetails > 0) {
        printf("largeShifts: ops=%i, sumNcol=%i, maxNcol=%i\n",
               eigenData.ops, eigenData.sumNcol, eigenData.maxNcol);
        printf("largeShifts: %i matrix-vector ops, "
               "%i reorthogonalizations in %.2f sec (%li wall)\n",
               info.matvec, info.north,
               (clock()-t1)/(float) CLOCKS_PER_SEC, tf-t2);
        printf("largeShifts:  %.2fmv cpu sec for ops\n",
               eigenData.tcpu_mv/(float) CLOCKS_PER_SEC);
        fflush(stdout);
    }

    if(info.stat < 0 || info.nec < info.ned) {
        fprintf(stderr, "trlan exit with stat=%i, "
                "finding %i of %i eigenpairs.\n",
               info.stat, info.nec, info.ned);
        exit(8);
    }

    /* We have calculated eigenvalues for hess^2.
       Convert to abs(eigenvalues) of hess itself. */
    for(k=0; k<*nvals; k++) {
        (*eval)[k] = sqrt((*eval)[k]);
    }
}

/* The extra parameter mvparam is not used in this case. */
void hessOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(columns(eigenData.matrix) == nrow);
    assert(rows(eigenData.matrix) == nrow);
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
void hessOp2(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(columns(eigenData.matrix) == nrow);
    assert(rows(eigenData.matrix) == nrow);
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
    int i;

    y = malloc(rows(hess) * sizeof(double));
    hessOp(rows(hess), 1, grad, 1, y, 1, NULL);
    printf("Dynamic part of hess.grad\n");
    for(i=0; i<rows(hess); i++) {
        printf("  %le\n", y[i]);
    }
    free(y);
}
