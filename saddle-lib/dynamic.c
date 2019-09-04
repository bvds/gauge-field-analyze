#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
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
} gaugeData;

void dynamicInit(unsigned int n, SparseRow *gauge,
		 unsigned int gaugeDimension, unsigned int gaugeElements,
		 cJSON *options) {
    gaugeData.n = n;
    gaugeData.matrix = gauge;
    gaugeData.matrixDimension = gaugeDimension;
    gaugeData.matrixElements = gaugeElements;
    gaugeData.options = options;
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

    if(istop >= 7) {
        printf("MINRES returned with istop=%i in %s, exiting.\n", istop, __FILE__);
        exit(7);
    }
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
    int k;
    memset(out, 0, n*sizeof(double));
    for(k=0; k<na; k++) {
        /* printf("Step:  %i %i %i %i %e %e %e\n", n, k, a[k].i, a[k].j, 
           a[k].value, in[a[k].j], out[a[k].i]);
           fflush(stdout); fflush(stderr); */
        out[a[k].i] += a[k].value * in[a[k].j];
    }
}

void vectorMatrix(const int n, const SparseRow *a, const int na,
                  const double *in, double *out) {
    int k;
    memset(out, 0, n*sizeof(double));
    for(k=0; k<na; k++) {
        out[a[k].j] += a[k].value * in[a[k].i];
    }
}
