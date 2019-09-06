/*
  Find eigenvalues associated with large shifts.
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

SparseMatrix *hessData;

void hessInit(SparseMatrix *hess) {
    hessData = hess;
}


void largeShifts(int n, double *initialVector, cJSON *options,
		 double **eval, double **evec, unsigned int *nvals) {
    // Let trlan figure out the work array allocation.
    const int lwrk = 0;
    const int iguess = 1; // Use supplies the intial vector
    double *wrk = NULL;
    int mev, maxlan, lohi, ned, maxmv;
    double tol = -1.0; // Use default tolerance: sqrt(machine epsilon)
    int nrow = n; // number of rows on this processor
    trl_info info;
    int i, printDetails = 1, restart = 1;
    cJSON *tmp;
    clock_t t1;
    time_t t2, tf;

    t1 = clock();
    time(&t2);

    // Would need a separate flag for the lohi == 0 case.
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenPairs");
    if(cJSON_IsNumber(tmp)) {
        ned = abs(tmp->valueint);
        /* In the source code, lohi=-2 sorts by distance from 
           info -> ref (which is set to zero). */
        lohi = tmp->valueint<0?-2:1;
    } else {
        ned = 1;
        lohi = -2;
    }
    /* We interpret maxIterations to be number of matrix-vector multiples.
       This sets an upper limit on time used. */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    if(cJSON_IsNumber(tmp)) {
        maxmv = tmp->valueint;
    } else {
        maxmv = n*ned;  // Documented efault (maxmv<0) seems to be broken.
    }
    /* maxLanczos is maximum number of Lanczos vectors to use.
       This sets an upper bound on the memory used. */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxLanczos");
    if(cJSON_IsNumber(tmp)) {
        maxlan = tmp->valueint;
    } else {
        maxlan = n;  // Since we do reorthogonalization.
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

    assert(n == hessData->columns);
    assert(n == hessData->rows);
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
    trlan(hessOp, dynamicProject, &info, nrow, mev, *eval, *evec, nrow, lwrk, wrk);
    *nvals = info.nec;
    if(printDetails > 1) {
        // This estimate of flops doesn't include dynamicProject() call. 
        trl_print_info(&info, hessData->nonzeros);
    } else if(printDetails > 0) {
        trl_terse_info(&info, stdout);
    }

    time(&tf);
    printf("largeShifts in %.2f sec (%li wall)\n",
           (clock()-t1)/(float) CLOCKS_PER_SEC, tf-t2);

    if(info.stat < 0 || info.nec < info.ned) {
        printf("trlan exit with stat=%i, finding %i of %i eigenpairs.\n",
               info.stat, info.nec, info.ned);
        exit(8);
    }
}

/* The extra parameter mvparam is not used in this case. */
void hessOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(hessData->columns == nrow);
    assert(hessData->rows == nrow);
    assert(mvparam == NULL);
    int k;

    for(k=0; k<ncol; k++) {
        matrixVector(hessData, xin+k*ldx, yout+k*ldy);
        // Use the fact that "in" and "out" can overlap.
        dynamicProject(nrow, yout+k*ldy, yout+k*ldy);
    }
}

// Debug print of dynamic part of hess.grad
void testOp(const int n, double *grad) {
    double *y;
    int i;
    y = malloc(n*sizeof(double));
    hessOp(n, 1, grad, 1, y, 1, NULL);
    printf("Dynamic part of hess.grad\n");
    for(i=0; i<n; i++) {
        printf("  %le\n", y[i]);
    }
    free(y);
}
