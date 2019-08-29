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
    unsigned int n;
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
void dynamicProject(unsigned int n, double *in, double *out) {
    unsigned int i;
    doublereal *b = malloc(gaugeData.matrixDimension*sizeof(double));
    doublereal *x = malloc(gaugeData.matrixDimension*sizeof(double));
    doublereal shift = 0.0, maxxnorm = 1.0e4;
    doublereal trancond, acondlim = 1.0e15;
    logical disable = 0;
    integer nout = 6, istop, itn, itnlim;
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond, rtol;

    /* Currently:
       "rTolerance": 1e-12,
       "printDetails": false
     */
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "maxIterations");
    if(cJSON_IsNumber(tmp)) {
        itnlim = tmp->valueint;
    } else {
        itnlim = 4*n;  // Assuming no reorthogonalization
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        rtol = tmp->valuedouble;
    } else {
        rtol = 1.6e-15;  // Try pointer to NULL or explicit Machine Epsilon?
    }

    assert(n == gaugeData.n);
    gaugeData.z = malloc(n*sizeof(doublereal));

    matrixVector(gaugeData.matrixDimension, gaugeData.matrix,
                 gaugeData.matrixElements, in, b);

    /* Verify that supplying NULL or ((void *)0) as argument is 
       interpreted by Fortran as not supplying an argument. */
  
    trancond = acondlim; // Always use MINRES
    MINRESQLP(
              &gaugeData.matrixDimension, gaugeProduct, b, &shift, NULL,
              &disable, &nout, &itnlim, &rtol, &maxxnorm, &trancond, &acondlim,
              x, &istop, &itn, &rnorm, &arnorm, &xnorm, &anorm, &acond);

    vectorMatrix(n, gaugeData.matrix,
                 gaugeData.matrixElements, x, gaugeData.z);

    for(i=0; i<n; i++) {
        out[i] = in[i] - gaugeData.z[i];
    }

    free(gaugeData.z); gaugeData.z = NULL;
    free(b);
    free(x);
}

int gaugeProduct(integer *vectorLength, doublereal *x, doublereal *y) {
    assert(*vectorLength == gaugeData.matrixDimension);
    vectorMatrix(gaugeData.n, gaugeData.matrix,
                 gaugeData.matrixElements,
                 x, gaugeData.z);
    matrixVector(*vectorLength, gaugeData.matrix, gaugeData.matrixElements,
                 gaugeData.z, y);
    return 0;
}

/* in and out must be distinct */
void matrixVector(int n, SparseRow *a, int na, const double *in, double *out) {
    int k;
    memset(out, 0, n*sizeof(double));
    for(k=0; k<na; k++) {
        /* printf("Step:  %i %i %i %i %e %e %e\n", n, k, a[k].i, a[k].j, 
           a[k].value, in[a[k].j], out[a[k].i]);
           fflush(stdout); fflush(stderr); */
        out[a[k].i] += a[k].value * in[a[k].j];
    }
}

void vectorMatrix(int n, SparseRow *a, int na, const double *in, double *out) {
    int k;
    memset(out, 0, n*sizeof(double));
    for(k=0; k<na; k++) {
        out[a[k].j] += a[k].value * in[a[k].i];
    }
}
