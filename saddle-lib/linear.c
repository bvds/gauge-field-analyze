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
    SparseMatrix *matrix;
    doublereal *vecs;
    integer nvecs;
    doublereal *z;
} hessData;

struct {
    doublereal *z;
    doublereal *vecs;
    integer nvecs;
    int count;
    integer maxVecs;
    integer pointer;
} minresOrtho;

void linearInit(SparseMatrix *hess, double *vecs, int nvecs) {
    hessData.matrix = hess;
    hessData.vecs = vecs;
    hessData.nvecs = nvecs;
}

/* "in" and "out" may overlap */
void linearSolve(integer n, double *b, cJSON *options, double *x) {
    doublereal shift = 0.0, *maxxnormp = NULL;
    // Explicit value for these, so we can force MINRES alogrithm.
    doublereal trancond, acondlim = 1.0e15;
    logical *disablep = NULL;
    integer nout, *noutp = NULL, istop, itn,
        itnlim = n; // Assuming reorthogonalization
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond,
        rtol, *rtolp = NULL, *abstolp = NULL;
    int printDetails = 0;
    clock_t t1;
    time_t t2, tf;

    t1 = clock();
    time(&t2);

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    if(cJSON_IsNumber(tmp)) {
        itnlim = tmp->valueint;
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxLanczosVecs");
    minresOrtho.maxVecs = cJSON_IsNumber(tmp)?tmp->valueint:itnlim;
    //
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        rtol = tmp->valuedouble;
        rtolp = &rtol;
    }
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "printDetails");
    if(cJSON_IsNumber(tmp)) {
        printDetails = tmp->valueint;
    } else if(cJSON_IsBool(tmp)) {
        printDetails = cJSON_IsTrue(tmp)?1:0;
    }
    if(printDetails>1) {
        nout = 6;  // Write to stdout
        noutp = &nout;
        printf("linearSolve writing to stdout\n");
    }

    assert(n>=0);
    assert(abs(n) == hessData.matrix->columns);
    assert(abs(n) == hessData.matrix->rows);
    hessData.z = malloc(hessData.nvecs*sizeof(doublereal));
    if((minresOrtho.z = malloc(minresOrtho.maxVecs*
                               sizeof(doublereal))) == NULL) {
        fprintf(stderr, "Can't allocated minresOrtho.z size %i*%li\n",
                minresOrtho.maxVecs, sizeof(doublereal));
        exit(124);
    }
    if((minresOrtho.vecs = malloc(minresOrtho.maxVecs*n*
                                  sizeof(doublereal))) == NULL) {
        fprintf(stderr, "Can't allocated minresOrtho.vecs size %i*%i*%li\n",
                minresOrtho.maxVecs, n, sizeof(doublereal));
        exit(125);
    }
    minresOrtho.nvecs = 0;
    minresOrtho.pointer = 0;
    minresOrtho.count = 0;

    /* grad may have components in the large shift direction */
    largeShiftProject(n, b);

    /* Verify that supplying NULL or ((void *)0) as argument is 
       interpreted by Fortran as not supplying an argument. */

    trancond = acondlim; // Always use MINRES
    MINRESQLP(
              &n, hessProduct, b, &shift, NULL,
	      (S_fp) userOrtho,
              disablep, noutp, &itnlim, rtolp, abstolp, maxxnormp,
              &trancond, &acondlim,
              x, &istop, &itn, &rnorm, &arnorm, &xnorm, &anorm, &acond);

    free(hessData.z);
    free(minresOrtho.z);
    free(minresOrtho.vecs);
    minresOrtho.nvecs = 0;

    time(&tf);
    if(printDetails > 0) {
        printf("linearSolve:  %i iterations, %i reorthogonalizations "
               "in %.2f sec (%li wall)\n",
               itn, minresOrtho.count,
               (clock()-t1)/(float) CLOCKS_PER_SEC, tf-t2);
        fflush(stdout);
    }

    if(istop >= 7) {
        fprintf(stderr, "MINRES returned with istop=%i in %s, exiting.\n", istop, __FILE__);
        exit(14);
    }
}

void hessProduct(integer *n, doublereal *x, doublereal *y) {
    assert(abs(*n) == hessData.matrix->rows);
    assert(abs(*n) == hessData.matrix->columns);
    hessOp(*n, 1, x, *n, y, *n, NULL);
    largeShiftProject(*n, y);
}

// Project out vecs, assuming they are orthonormal
void userOrtho(char *action, integer *n, double *y) {
    const char trans='T', normal='N';
    const double one=1.0, zero=0.0, minusone=-1.0;
    const integer inc=1;
    const integer dn = *n * sizeof(double);

    assert(abs(*n) == hessData.matrix->rows);
    assert(abs(*n) == hessData.matrix->columns);

    if(*action=='a') {
        /* add vector to ortho list */
        memcpy(minresOrtho.vecs + minresOrtho.pointer*(*n), y, dn);
        minresOrtho.pointer += 1;
        if(minresOrtho.pointer==minresOrtho.maxVecs) {
            minresOrtho.pointer = 0;
        }
        if(minresOrtho.nvecs<minresOrtho.maxVecs) {
            minresOrtho.nvecs += 1;
        }
    } else if(*action=='o') {
        /* orthogonalize vector against previous
           vectors, gauge transforms, and large shifts 

           Need to figure out the best order for
           the three orthogonalizations.
        */
        minresOrtho.count += 1;

        DGEMV(&trans, n, &minresOrtho.nvecs, &one,
              minresOrtho.vecs, n, y, &inc, &zero,
              minresOrtho.z, &inc);
        DGEMV(&normal, n, &minresOrtho.nvecs, &minusone,
              minresOrtho.vecs, n, minresOrtho.z, &inc, &one,
              y, &inc);

        largeShiftProject(*n, y);
        
        dynamicProject(*n, y, NULL);
    } else {
        fprintf(stderr, "Bad value %s\n", action);
        exit(111);
    }
}

void largeShiftProject(integer n, double *y) {
    const char trans='T', normal='N';
    const double one=1.0, zero=0.0, minusone=-1.0;
    const integer inc=1;

    DGEMV(&trans, &n, &hessData.nvecs, &one,
          hessData.vecs, &n, y, &inc, &zero,
          hessData.z, &inc);
    DGEMV(&normal, &n, &hessData.nvecs, &minusone,
          hessData.vecs, &n, hessData.z, &inc, &one,
          y, &inc);
}
