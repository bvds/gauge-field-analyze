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
    mat_int nrow;
    SparseMatrix *matrix;
    doublereal *vecs;
    integer nvecs;
    doublereal *z;
#ifdef USE_MPI
    MPI_Comm mpicom;
#endif
} hessData;

struct {
    doublereal *z;
    doublereal *vecs;
    integer nvecs;
    int count;
    integer maxVecs;
    integer pointer;
#ifdef USE_MPI
    MPI_Comm mpicom;
#endif
} minresOrtho;

void linearInit(SparseMatrix *hess, const mat_int nrow,
                double *vecs, int nvecs, void *mpicomp) {
    hessData.nrow = nrow;  // For sanity checks
    hessData.matrix = hess;
    hessData.vecs = vecs;
    hessData.nvecs = nvecs;
#ifdef USE_MPI
    hessData.mpicom = *((MPI_Comm *) mpicomp);
#else
    assert(mpicomp == NULL);
#endif
}

/* "in" and "out" may overlap */
void linearSolve(integer n, double *b, cJSON *options, double *x) {
    doublereal shift = 0.0, *maxxnormp = NULL;
    // Explicit value for these, so we can force MINRES alogrithm.
    doublereal trancond, acondlim = 1.0e15;
    logical *disablep = NULL;
    integer nout, *noutp = NULL, istop, itn, itnlim;
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond,
        rtol, *rtolp = NULL, *abstolp = NULL;
    int printDetails = 0, wrank;
    unsigned long tcpu;
    clock_t t1;
    time_t t2, tf;

    t1 = clock();
    time(&t2);

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    itnlim = cJSON_IsNumber(tmp)?(mat_int) tmp->valueint:
        // Default appropriate for full reorthogonalization
        hessData.matrix->rows;
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

    assert(hessData.nrow == (mat_int) n);
    hessData.z = malloc(hessData.nvecs*sizeof(doublereal));
    if((minresOrtho.z = malloc(minresOrtho.maxVecs*
                               sizeof(doublereal))) == NULL) {
        fprintf(stderr, "Can't allocate minresOrtho.z size %i*%li\n",
                minresOrtho.maxVecs, sizeof(doublereal));
        exit(124);
    }
    if((minresOrtho.vecs = malloc(minresOrtho.maxVecs*n*
                                  sizeof(doublereal))) == NULL) {
        fprintf(stderr, "Can't allocate minresOrtho.vecs size %i*%i*%li\n",
                minresOrtho.maxVecs, n, sizeof(doublereal));
        exit(125);
    }
    minresOrtho.nvecs = 0;
    minresOrtho.pointer = 0;
    minresOrtho.count = 0;
#ifdef USE_MPI
    minresOrtho.mpicom = hessData.mpicom;
#endif

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
    tcpu = clock() - t1;
#ifdef USE_MPI
    MPI_Comm_rank(hessData.mpicom, &wrank);
    MPI_Allreduce(MPI_IN_PLACE, &tcpu, 1, MPI_LONG,
                  MPI_SUM, hessData.mpicom);
#else
    wrank = 0;
#endif
    if(printDetails > 0 && wrank==0) {
        printf("linearSolve:  %i iterations, %i reorthogonalizations "
               "in %.2f sec (%li wall)\n",
               itn, minresOrtho.count,
               tcpu/(float) CLOCKS_PER_SEC, tf-t2);
        fflush(stdout);
    }

    if(istop >= 7) {
        fprintf(stderr, "%i:  MINRES returned with istop=%i in %s, exiting.\n",
                wrank, istop, __FILE__);
        exit(14);
    }
}

void hessProduct(integer *n, doublereal *x, doublereal *y) {
    hessOp(*n, 1, x, *n, y, *n, NULL);
    largeShiftProject(*n, y);
}

// Project out vecs, assuming they are orthonormal
void userOrtho(char *action, integer *n, double *y) {
    const char trans='T', normal='N';
    const double one=1.0, zero=0.0, minusone=-1.0;
    const integer inc=1;
    const integer dn = *n * sizeof(double);

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
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, minresOrtho.z, minresOrtho.nvecs,
                      MPI_DOUBLE_PRECISION, MPI_SUM, minresOrtho.mpicom);
#endif
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

    assert(hessData.nrow == (mat_int) n);

    DGEMV(&trans, &n, &hessData.nvecs, &one,
          hessData.vecs, &n, y, &inc, &zero,
          hessData.z, &inc);
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, hessData.z, hessData.nvecs,
                  MPI_DOUBLE_PRECISION, MPI_SUM, hessData.mpicom);
#endif
    DGEMV(&normal, &n, &hessData.nvecs, &minusone,
          hessData.vecs, &n, hessData.z, &inc, &one,
          y, &inc);
}
