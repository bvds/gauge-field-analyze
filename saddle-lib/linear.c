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

struct {
    doublereal *z;
    doublereal *vecs;
    integer nvecs;    
} minresOrtho;

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
    minresOrtho.z = malloc(itnlim*sizeof(doublereal));
    minresOrtho.vecs = NULL; minresOrtho.nvecs = 0;

    /* grad may have components in the large shift direction */
    largeShiftProject(n, b);

    /* Verify that supplying NULL or ((void *)0) as argument is 
       interpreted by Fortran as not supplying an argument. */

    trancond = acondlim; // Always use MINRES
    MINRESQLP(
              &n, hessProduct, b, &shift, NULL, (S_fp) userOrtho,
              &disable, &nout, &itnlim, &rtol, &maxxnorm, &trancond, &acondlim,
              x, &istop, &itn, &rnorm, &arnorm, &xnorm, &anorm, &acond);

    free(hessData.z);
    free(minresOrtho.z);
    minresOrtho.vecs = realloc(minresOrtho.vecs, 0);
    minresOrtho.nvecs = 0;
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
void userOrtho(char *action, integer *n, double *y) {
    const char trans='T', normal='N';
    const double one=1.0, zero=0.0, minusone=-1.0;
    const integer inc=1;
    integer dn = *n * sizeof(double);

    printf("userOrtho n=%i hessData.n=%i\n", *n, hessData.n);
    assert(*n == hessData.n);

    if(*action=='a') {
        /* add vector to ortho list */
        minresOrtho.vecs = realloc(minresOrtho.vecs, (minresOrtho.nvecs+1) * dn);
        minresOrtho.z = realloc(minresOrtho.z, (minresOrtho.nvecs + 1) * sizeof(double));
        memcpy(y, minresOrtho.vecs + minresOrtho.nvecs * dn, dn);
        minresOrtho.nvecs += 1;
    } else if(*action=='o') {
        /* orthogonalize vector against previous
           vectors, gauge transforms, and large shifts 

           Need to figure out the best order for
           the three orthogonalizations.
        */

        // BvdS:  dimensions/transposes are all screwed up.
        DGEMV(&trans, n, &minresOrtho.nvecs, &one,
              minresOrtho.vecs, n, y, &inc, &zero,
              minresOrtho.z, &inc);
        DGEMV(&normal, n, &minresOrtho.nvecs, &minusone,
              minresOrtho.vecs, n, minresOrtho.z, &inc, &one,
              y, &inc);

        largeShiftProject(*n, y);
        
        dynamicProject(*n, y, y);
    } else {
        printf("Bad value %s\n", action);
        exit(1);
    }
}

void largeShiftProject(integer n, double *y) {
    const char trans='T', normal='N';
    const double one=1.0, zero=0.0, minusone=-1.0;
    const integer inc=1;

    // BvdS:  dimensions/transposes are all screwed up.
    DGEMV(&trans, &n, &hessData.nvecs, &one,
          hessData.vecs, &n, y, &inc, &zero,
          hessData.z, &inc);
#if 0
    int i;
    for(i=0; i<hessData.nvecs; i++){
        printf("z(%i)=%lf\n", i, hessData.z[i]);
    }
#endif
    DGEMV(&normal, &n, &hessData.nvecs, &minusone,
          hessData.vecs, &n, hessData.z, &inc, &one,
          y, &inc);
}
