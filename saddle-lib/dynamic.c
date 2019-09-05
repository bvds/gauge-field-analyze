#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "shifts.h"
/*
  Project out infinitesimal gauge transforms from
  the shifts using MINRES/MINRES-QLP
*/

struct {
    int n;
    SparseRow *matrix;
    integer matrixDimension;
    unsigned int matrixElements;
    cJSON *options;
    doublereal *z;
    unsigned long int tcpu;
    unsigned long int twall;
    unsigned int count;
    unsigned int matVec;
} gaugeData;

void dynamicInit(unsigned int n, SparseRow *gauge,
		 unsigned int gaugeDimension, unsigned int gaugeElements,
		 cJSON *options) {
    gaugeData.n = n;
    gaugeData.matrix = gauge;
    gaugeData.matrixDimension = gaugeDimension;
    gaugeData.matrixElements = gaugeElements;
    gaugeData.options = options;
    gaugeData.tcpu = 0;
    gaugeData.twall = 0;
    gaugeData.count = 0;
    gaugeData.matVec = 0;
}

/* "in" and "out" may overlap */
void dynamicProject(const int n, double *in, double *out) {
    int i;
    doublereal *b = malloc(gaugeData.matrixDimension*sizeof(double));
    doublereal *x = malloc(gaugeData.matrixDimension*sizeof(double));
    doublereal shift = 0.0, *maxxnormp = NULL;
    // Explicit value for these, so we can force MINRES alogrithm.
    doublereal trancond, acondlim = 1.0e15;
    logical *disablep = NULL;
    integer nout, *noutp = NULL, istop, itn, itnlim;
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond,
        rtol, *rtolp = NULL;
    clock_t t1;
    time_t t2, tf;

    t1 = clock();
    time(&t2);

    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "maxIterations");
    if(cJSON_IsNumber(tmp)){
        itnlim = tmp->valueint;
    } else {
        // Appropriate for reorthogonalization
        itnlim = n;
    }
    /* 
       We may want to adjust rtol in cases where norm(b) is very small,
       since this is often just removing roundoff error.

       In any case, the gauge configurations are single precision,
       so rtol is normally about 1e-7.
    */
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        rtol = tmp->valuedouble;
        rtolp = &rtol;
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "printDetails");
    if((cJSON_IsNumber(tmp) && tmp->valueint>0) || cJSON_IsTrue(tmp)) {
        nout = 6;  // Write to stdout
        noutp = &nout;
    }

    assert(n == gaugeData.n);
    gaugeData.z = malloc(n*sizeof(doublereal));

    matrixVector(gaugeData.matrixDimension, gaugeData.matrix,
                 gaugeData.matrixElements, in, b);
  
    trancond = acondlim; // Always use MINRES
    MINRESQLP(&gaugeData.matrixDimension, gaugeProduct, b, &shift, NULL, NULL,
              disablep, noutp, &itnlim, rtolp, maxxnormp, &trancond, &acondlim,
              x, &istop, &itn, &rnorm, &arnorm, &xnorm, &anorm, &acond);

    vectorMatrix(n, gaugeData.matrix,
                 gaugeData.matrixElements, x, gaugeData.z);

    for(i=0; i<n; i++) {
        out[i] = in[i] - gaugeData.z[i];
    }

    free(gaugeData.z); gaugeData.z = NULL;
    free(b);
    free(x);

    time(&tf);
    gaugeData.tcpu += clock()-t1;
    gaugeData.twall += tf - t2; 
    gaugeData.count += 1;
    gaugeData.matVec += itn;

    if(istop >= 7) {
        printf("MINRES returned with istop=%i in %s, exiting.\n", istop, __FILE__);
        exit(7);
    }
}

void printDynamicStats() {
    printf("dynamicProject %i calls (%i iterations) in %.2f sec (%li wall)\n",
           gaugeData.count, gaugeData.matVec,
           gaugeData.tcpu/(float) CLOCKS_PER_SEC, gaugeData.twall);
}

int gaugeProduct(const integer *vectorLength, const doublereal *x,
                 doublereal *y) {
    assert(*vectorLength == gaugeData.matrixDimension);
    vectorMatrix(gaugeData.n, gaugeData.matrix,
                 gaugeData.matrixElements,
                 x, gaugeData.z);
    matrixVector(*vectorLength, gaugeData.matrix, gaugeData.matrixElements,
                 gaugeData.z, y);
    return 0;
}

/* in and out must be distinct */
void matrixVector(const int n, const SparseRow *a, const int na,
                  const double *in, double *out) {
#ifdef USE_LIBRSB
    const rsb_trans_t trans = RSB_TRANSPOSITION_N;
    const double zero = 0.0, one = 1.0;
    int errval;

    // Sanity test
    int nna;
    rsb_mtx_get_info(a, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &nna);
    assert(n == nna);
    rsb_mtx_get_info(a, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &nna);
    // Only one triangle of a symmetric matrix is stored.
    assert(na == 2*nna-n || na == nna);

    if((errval = rsb_spmv(trans, &one, a, in, 1, &zero, out, 1))
       != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_spmv error 0x%x, exiting\n", errval);
        exit(121);
    }
#if 0  // Debug print
    double norm;
    int k;
    if((errval = rsb_mtx_get_nrm(a, &norm, RSB_EXTF_NORM_TWO))
       != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_mtx_get_nrm error 0x%x, exiting\n", errval);
        exit(121);
    }
    printf("== matrixVector for %i %i norm=%lf\n", n, nna, norm);
    for(k=0; k<n; k++) {
        printf("==   %le %le\n", in[k], out[k]);
    }
#endif
#else
    int k;
    memset(out, 0, n*sizeof(double));
    for(k=0; k<na; k++) {
        out[a[k].i] += a[k].value * in[a[k].j];
    }
#endif
}

void vectorMatrix(const int n, const SparseRow *a, const int na,
                  const double *in, double *out) {
#ifdef USE_LIBRSB
    const rsb_trans_t trans = RSB_TRANSPOSITION_T;
    const double zero = 0.0, one = 1.0;

    // Sanity test
    int nna;
    rsb_mtx_get_info(a, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T , &nna);
    assert(n == nna);
    rsb_mtx_get_info(a, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &nna);
    assert(na == nna);

    rsb_spmv(trans, &one, a, in, 1, &zero, out, 1);
#else
    int k;
    memset(out, 0, n*sizeof(double));
    for(k=0; k<na; k++) {
        out[a[k].j] += a[k].value * in[a[k].i];
    }
#endif
}
