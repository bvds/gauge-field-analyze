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

struct {
    SparseMatrix *matrix;
} eigenData;

void largeShifts(SparseMatrix *hess, double *initialVector, cJSON *options,
		 double **eval, double **evec, unsigned int *nvals) {
    // Let trlan figure out the work array allocation.
    const int lwrk = 0;
    const int iguess = 1; // Use supplies the intial vector
    double *wrk = NULL;
    int mev, maxlan, lohi, ned, maxmv;
    double tol = -1.0; // Use default tolerance: sqrt(machine epsilon)
    int nrow = rows(hess); // number of rows on this processor
    trl_info info;
    int i, printDetails = 1, restart = 1;
    cJSON *tmp;
    clock_t t1;
    time_t t2, tf;

    t1 = clock();
    time(&t2);

    // Used by HessOp
    eigenData.matrix = hess;

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
    trlan(hessOp, dynamicProject, &info, nrow, mev, *eval, *evec,
          nrow, lwrk, wrk);
    *nvals = info.nec;
    if(printDetails > 1) {
        /* This estimate of flops doesn't include dynamicProject() call. 
           Assuming 1 addition and 1 multiply per nonzero element. */
        trl_print_info(&info, 2*nonzeros(hess));
    } else if(printDetails > 0) {
        trl_terse_info(&info, stdout);
    }

    time(&tf);
    if(printDetails > 0) {
        printf("largeShifts: %i matrix-vector ops, "
               "%i reorthogonalizations in %.2f sec (%li wall)\n",
               info.matvec, info.north,
               (clock()-t1)/(float) CLOCKS_PER_SEC, tf-t2);
    }

    if(info.stat < 0 || info.nec < info.ned) {
        printf("trlan exit with stat=%i, finding %i of %i eigenpairs.\n",
               info.stat, info.nec, info.ned);
        exit(8);
    }
}

/* The extra parameter mvparam is not used in this case. */
void hessOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(columns(eigenData.matrix) == nrow);
    assert(rows(eigenData.matrix) == nrow);
    assert(mvparam == NULL);
    int k;

    for(k=0; k<ncol; k++) {
        matrixVector(eigenData.matrix, xin+k*ldx, yout+k*ldy);
        // Use the fact that "in" and "out" can overlap.
        dynamicProject(nrow, yout+k*ldy, yout+k*ldy);
    }
}

// Debug print of dynamic part of hess.grad
void testOp(SparseMatrix *hess, double *grad) {
    double *y;
    int i;

    eigenData.matrix = hess;

    y = malloc(rows(hess) * sizeof(double));
    hessOp(rows(hess), 1, grad, 1, y, 1, NULL);
    printf("Dynamic part of hess.grad\n");
    for(i=0; i<rows(hess); i++) {
        printf("  %le\n", y[i]);
    }
    free(y);
}
