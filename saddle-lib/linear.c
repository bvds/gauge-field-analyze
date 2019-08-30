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
    integer n;
    SparseRow *matrix;
    unsigned int matrixElements;
    doublereal *vecs;
    integer nvecs;
    doublereal *z;
} hessData;

void linearInit(unsigned int n, SparseRow *hess, int hessElements,
                double *vecs, int nvecs) {
    hessData.n = n;
    hessData.matrix = hess;
    hessData.matrixElements = hessElements;
    hessData.vecs = vecs;
    hessData.nvecs = nvecs;
}

/* "in" and "out" may overlap */
void linearSolve(integer n, double *b, cJSON *options, double *x) {
    doublereal shift = 0.0, maxxnorm = 1.0e4;
    doublereal trancond, acondlim = 1.0e15;
    logical disable = 0;
    integer nout = 6, istop, itn, itnlim;
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond, rtol;

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    if(cJSON_IsNumber(tmp)) {
        itnlim = tmp->valueint;
    } else {
        itnlim = 4*n;  // Assuming no reorthogonalization
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        rtol = tmp->valuedouble;
    } else {
        rtol = 1.6e-15;  // Try pointer to NULL or explicit Machine Epsilon?
    }

    assert(n == hessData.n);
    hessData.z = malloc(hessData.nvecs*sizeof(doublereal));

    /* grad may have components in the large shift direction */
    largeShiftProject(n, b);

    /* Verify that supplying NULL or ((void *)0) as argument is 
       interpreted by Fortran as not supplying an argument. */

    trancond = acondlim; // Always use MINRES
    MINRESQLP(
              &n, hessProduct, b, &shift, NULL,
              &disable, &nout, &itnlim, &rtol, &maxxnorm, &trancond, &acondlim,
              x, &istop, &itn, &rnorm, &arnorm, &xnorm, &anorm, &acond);

    free(hessData.z); hessData.z = NULL;
}

int hessProduct(integer *n, doublereal *x, doublereal *y) {
    assert(*n == hessData.n);
    matrixVector(*n, hessData.matrix, hessData.matrixElements, x, y);
    // Use the fact that "in" and "out" can overlap.
    dynamicProject(*n, y, y);
    largeShiftProject(*n, y);
    return 0;
}

// Project out vecs, assuming they are orthonormal
// "in" and "out" may overlap
void largeShiftProject(integer n, double *y) {
    const char yes='Y', no='N';
    const double one=1.0, zero=0.0, minusone=-1.0;
    const integer inc=1;

    assert(n == hessData.n);
    
    // BvdS:  dimensions/transposes are all screwed up.
    DGEMV(&yes, &hessData.nvecs, &n, &one,
          hessData.vecs, &n, y, &inc, &zero,
          hessData.z, &inc);
    DGEMV(&no, &hessData.nvecs, &n, &minusone,
          hessData.vecs, &n, hessData.z, &inc, &one,
          y, &inc);
}
