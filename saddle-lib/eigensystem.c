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
// #include "trl_map.h"
#include "trl_comm_i.h"

struct {
    int n;
    SparseRow *matrix;
    unsigned int matrixElements;
} hessData;

void hessInit(unsigned int n, SparseRow *hess, unsigned int hessElements) {
    hessData.n = n;
    hessData.matrix = hess;
    hessData.matrixElements = hessElements;
}


/* The extra parameter mvparam is not used in this case. */
void hessOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(hessData.n == nrow);
    assert(mvparam == NULL);
    int k;

    for(k=0; k<ncol; k++) {
        matrixVector(nrow, hessData.matrix,
                     hessData.matrixElements, xin+k*ldx, yout+k*ldy);
        // Use the fact that "in" and "out" can overlap.
        dynamicProject(nrow, yout+k*ldy, yout+k*ldy);
    }
}

void largeShifts(int n, double *initialVector, cJSON *options,
		 double **eval, double **evec, unsigned int *nvals) {
    int lwrk, mev, maxlan, lohi, ned, maxmv;
    double *wrk, rtol;
    int nrow = n; // number of rows on this processor
    trl_info info;
    int i;
    cJSON *tmp;

    /*
       Currently defined:
        "eigenPairs": -35,
        "maxIterations": 3000,
        "innerInterval": 10,
        "rTolerance": 1e-7,
        "printDetails": 2,
        "printsPerDecade": 10,
        "restartLanczos": false
    */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenPairs");
    if(cJSON_IsNumber(tmp)) {
        mev = abs(tmp->valueint);
        lohi = tmp->valueint;
    } else {
        mev = 1;
        lohi = -1;
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    if(cJSON_IsNumber(tmp)) {
        maxlan = tmp->valueint;
    } else {
        maxlan = 4*n;  // Assuming no reorthogonalization
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        rtol = tmp->valuedouble;
    } else {
        rtol = 1.6e-8;  // Should try NULL or explicit sqrt(Machine Epsilon)
    }

    assert(n == hessData.n);

    *eval = malloc(mev*sizeof(double));
    *evec = malloc(mev*nrow*sizeof(double));
    // initialize info -- tell TRLAN to compute NEV smallest eigenvalues
    // of the diag_op
    lwrk=maxlan*(maxlan+10);
    if( lwrk > 0 ) {
	wrk = malloc(lwrk*sizeof(double));
    }
    ned = mev; // BvdS:  not sure what the difference is
    maxmv = ned*nrow;  // BvdS:  max number of matrix multiplies.
    trl_init_info(&info, nrow, maxlan, lohi, ned, rtol, 1, maxmv, 0);
    trl_set_iguess(&info, 0, 1, 0, NULL);
    // The Lanczos recurrence is set to start with initial guess
    memset(eval, 0, mev*sizeof(double));
    for( i=0; i<nrow; i++ ) {
        *evec[i] = initialVector==NULL?1.0:initialVector[i];
    }
    //info.verbose =  8;
    // call TRLAN to compute the eigenvalues
    trlan(hessOp, &info, nrow, mev, *eval, *evec, nrow, lwrk, wrk);
    trl_print_info(&info, 3*nrow);
    *nvals = 0;

    if( lwrk > 0 ) {
	free(wrk);
    }
}
