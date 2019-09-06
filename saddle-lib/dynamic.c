#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "shifts.h"
/*
  Project out infinitesimal gauge transforms from
  the shifts using MINRES/MINRES-QLP
*/

struct {
    SparseMatrix *matrix;
    cJSON *options;
    doublereal *z;
    unsigned long int tcpu;
    unsigned long int twall;
    unsigned int count;
    unsigned int matVec;
} gaugeData;

void dynamicInit(SparseMatrix *gauge, cJSON *options) {
    gaugeData.matrix = gauge;
    gaugeData.options = options;
    gaugeData.tcpu = 0;
    gaugeData.twall = 0;
    gaugeData.count = 0;
    gaugeData.matVec = 0;

    /* Calculate the largest and smallest eigenvalues
       of gaugeProduct.  These will inform the stopping condition
       for MINRES. */
    
}

/* "in" and "out" may overlap */
void dynamicProject(const int n, double *in, double *out) {
    int i;
    doublereal *b = malloc(gaugeData.matrix->rows * sizeof(double));
    doublereal *x = malloc(gaugeData.matrix->rows * sizeof(double));
    doublereal shift = 0.0, *maxxnormp = NULL;
    // Explicit value for these, so we can force MINRES alogrithm.
    doublereal trancond, acondlim = 1.0e15;
    logical *disablep = NULL;
    integer nout, *noutp = NULL, istop, itn, itnlim;
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond,
        rtol, *rtolp = NULL, *abstolp = NULL;
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

    assert(n == gaugeData.matrix -> columns);
    gaugeData.z = malloc(n*sizeof(doublereal));

    matrixVector(gaugeData.matrix, in, b);
  
    trancond = acondlim; // Always use MINRES
    MINRESQLP(&(gaugeData.matrix->rows), gaugeProduct, b, &shift, NULL, NULL,
              disablep, noutp, &itnlim, rtolp, abstolp, maxxnormp, &trancond, &acondlim,
              x, &istop, &itn, &rnorm, &arnorm, &xnorm, &anorm, &acond);

    vectorMatrix(gaugeData.matrix, x, gaugeData.z);
    // This could be combined with the matrix product, BLAS style.
    for(i=0; i<n; i++) {
        out[i] = in[i] - gaugeData.z[i];
    }

    double norm2x = 0.0, norm2z = 0.0, norm2b = 0.0; //, indotout = 0.0; 
    for(i=0; i<n; i++) {
        norm2x += out[i]*out[i];
        norm2z += gaugeData.z[i]*gaugeData.z[i];
        // indotout += out[i]*gaugeData.z[i];
    }
    for(i=0; i<gaugeData.matrix->rows; i++) {
        norm2b += b[i]*b[i];
    }
#if 0
    printf("** dynamicProject norm(x)=%le shift ratio=%le norm(b)=%le itn=%i\n",
           sqrt(norm2x), sqrt(norm2z/norm2x), sqrt(norm2b), itn);
           //indotout/sqrt(norm2x*norm2z));
#endif

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
    assert(*vectorLength == gaugeData.matrix->rows);
    vectorMatrix(gaugeData.matrix, x, gaugeData.z);
    matrixVector(gaugeData.matrix, gaugeData.z, y);
    return 0;
}

/* in and out must be distinct */
void matrixVector(const SparseMatrix *a,
                  const double *in, double *out) {
#ifdef USE_LIBRSB
    const rsb_trans_t trans = RSB_TRANSPOSITION_N;
    const double zero = 0.0, one = 1.0;
    int errval;

    if((errval = rsb_spmv(trans, &one, a, in, 1, &zero, out, 1))
       != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_spmv error 0x%x, exiting\n", errval);
        exit(121);
    }
#else
    unsigned int k;
    SparseRow *row;
    memset(out, 0, a->rows * sizeof(double));
    for(k=0; k<a->nonzeros; k++) {
        row = a->data+k;
        out[row->i] += row->value * in[row->j];
    }
#endif
}

void vectorMatrix(const SparseMatrix *a,
                  const double *in, double *out) {
#ifdef USE_LIBRSB
    const rsb_trans_t trans = RSB_TRANSPOSITION_T;
    const double zero = 0.0, one = 1.0;
    int errval;

    if((errval = rsb_spmv(trans, &one, a, in, 1, &zero, out, 1))
       != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_spmv error 0x%x, exiting\n", errval);
        exit(121);
    }
#else
    unsigned int k;
    SparseRow *row;
    memset(out, 0, a->columns * sizeof(double));
    for(k=0; k<a->nonzeros; k++) {
        row = a->data+k;
        out[row->j] += row->value * in[row->i];
    }
#endif
}
