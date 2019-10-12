#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "shifts.h"
#include "trlan.h"
#include "trl_comm_i.h"

/*
  Project out infinitesimal gauge transforms from
  the shifts using MINRES/MINRES-QLP applied to the
  square of the infinitesimal gauge transform shift matrix:

      v <- v - o^T.x

   where

      A.x = b, A = o.o^T, b = o.x

  In principle, one could apply a linear least squares
  solver to o itself.  However, A is pretty well-conditioned
  so the conjugate gradient solution is relatively fast.

  Also, we terminate the Conjugate Gradient solver based
  the upper limit of norm(o^T.x) relative to norm(v).
*/

struct {
    SparseMatrix *matrix;
    cJSON *options;
    doublereal *z;
    doublereal *x;
    doublereal *b;
    unsigned long int tcpu;
    unsigned long int tcpu_mv;
    unsigned long int tcpu_vm;
    unsigned long int twall;
    mat_int count;
    mat_int usertol;
    mat_int matVec;
    int maxItn;
    int printDetails;
    doublereal minEigenvalue;
    // doublereal maxEigenvalue;
} gaugeData;



void dynamicInit(SparseMatrix *gauge, cJSON *options, void *mpicomp) {
    // Let trlan figure out the work array allocation.
    const int lwrk = 0; double *wrk = NULL;
    int mev, maxlan, lohi, ned, maxmv;
    double tol = -1.0; // Use default tolerance: sqrt(machine epsilon)
    double *eval, *evec;
    int nrow; // number of rows on this processor
    int restart = 3; // restart strategy
    trl_info info;
    cJSON *tmp;
    int i;

    gaugeData.matrix = gauge;
    gaugeData.options = options;
    gaugeData.tcpu = 0;
    gaugeData.twall = 0;
    gaugeData.count = 0;
    gaugeData.usertol = 0;
    gaugeData.matVec = 0;
    gaugeData.maxItn = 0;
    gaugeData.printDetails = 0;

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "printDetails");
    if(cJSON_IsNumber(tmp)) {
        gaugeData.printDetails = tmp->valueint;
    } else if(cJSON_IsBool(tmp)) {
        gaugeData.printDetails = cJSON_IsTrue(tmp)?1:0;
    }

#ifdef USE_MKL
    gaugeData.b = mkl_malloc(gauge->rows*sizeof(double), MALLOC_ALIGN);
    gaugeData.x = mkl_malloc(gauge->rows*sizeof(double), MALLOC_ALIGN);
    gaugeData.z = mkl_malloc(gauge->columns*sizeof(*gaugeData.z), MALLOC_ALIGN);
#else
    gaugeData.b = malloc(gauge->rows*sizeof(double));
    gaugeData.x = malloc(gauge->rows*sizeof(double));
    gaugeData.z = malloc(gauge->columns*sizeof(*gaugeData.z));
#endif

    /* Calculate the smallest eigenvalue of gaugeProduct.  
       This will inform the stopping condition for MINRES.

       Tried finding both the smallest and largest eigenvalues,
       but had very slow convergence for large matrices.
    */
    lohi = -1;
    nrow = gauge->rows;
    ned = 1;
    // Uses same maxIterations as MINRES.
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "maxIterations");
    maxmv = cJSON_IsNumber(tmp)?tmp->valueint:ned*nrow;
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "maxLanczosVecs");
    maxlan = cJSON_IsNumber(tmp)?tmp->valueint:nrow;

    mev = ned; // Allocate memory for the number of requested eigenpairs
    eval = malloc(mev*sizeof(double));
    evec = malloc(mev*nrow*sizeof(double));
    trl_init_info(&info, nrow, maxlan, lohi, ned, tol, restart, maxmv, mpicomp);
    memset(eval, 0, mev*sizeof(double));

    // call TRLAN to compute the eigenvalues
    trlan(gaugeOp, NULL, &info, nrow, mev, eval, evec, nrow, lwrk, wrk);
    if(gaugeData.printDetails > 1) {
        // 2 matrix multplies, with 1 add and 1 multiply per nonzero element.
        trl_print_info(&info, 4*(gauge->nonzeros));
    } else if(gaugeData.printDetails > 0) {
        trl_terse_info(&info, stdout);
    }
    if(gaugeData.printDetails > 0) {
        printf("Gauge matrix extreme eigenvalues found:\n");
        for(i=0; i<info.nec; i++) {
            printf("  %le\n", eval[i]);
        }
    }

    if(info.nec>0 && info.stat==0) {
        gaugeData.minEigenvalue = eval[0];
        // gaugeData.maxEigenvalue = eval[info.nec-1];
    } else {
        printf("trlan exit with stat=%i, finding %i of %i eigenpairs.\n",
               info.stat, info.nec, info.ned);
        exit(8);
    }

    free(eval); free(evec);
}

void dynamicProject(const integer n, double *v, double *normDiff) {
    doublereal shift = 0.0, *maxxnormp = NULL;
    // Explicit value for these, so we can force MINRES alogrithm.
    doublereal trancond, acondlim = 1.0e15;
    logical *disablep = NULL;
    integer rr, nout, *noutp = NULL, istop, itn, itnlim;
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond,
        rtol, *rtolp = NULL, abstol, *abstolp = NULL;
    doublereal normv;
    const integer one=1;
    const doublereal minusone = -1.0;
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
       The gauge configurations are single precision,
       so rtol is normally about 1e-7.
    */
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        rtol = tmp->valuedouble;
        rtolp = &rtol;
    }
    /*
      Tolerance of the projection relative to the initial vector.
      For now, tie this to rTolerance.  May want separate flag?
    */
    if(rtolp != NULL) {
        normv = DNRM2(&n, v, &one);
        abstol = rtol*normv*sqrt(gaugeData.minEigenvalue);
        abstolp = &abstol;
#if 0
        printf("Setting abstolp=%le normin=%le val=%le sqrt=%le\n",
               abstol, normin, gaugeData.minEigenvalue,
               sqrt(gaugeData.minEigenvalue));
#endif
    }
    if(gaugeData.printDetails > 2) {
        nout = 6;  // Write to stdout
        noutp = &nout;
    }

    // Sanity test for gaugeData.z
    assert(abs(n) == gaugeData.matrix->columns);

    matrixVector(gaugeData.matrix, v, gaugeData.b);

    trancond = acondlim; // Always use MINRES
    rr = gaugeData.matrix->rows;
    MINRESQLP(&rr, gaugeProduct, gaugeData.b,
              &shift, NULL, NULL, disablep, noutp, &itnlim,
              rtolp, abstolp, maxxnormp, &trancond, &acondlim,
              gaugeData.x, &istop, &itn, &rnorm, &arnorm,
              &xnorm, &anorm, &acond);

    vectorMatrix(gaugeData.matrix, gaugeData.x, gaugeData.z);
    // This could be combined with the above matrix product, BLAS style.
    DAXPY(&n, &minusone, gaugeData.z, &one, v, &one);
    if(normDiff != NULL){
        *normDiff = DNRM2(&n, gaugeData.z, &one);
    }

    time(&tf);
    gaugeData.tcpu += clock()-t1;
    gaugeData.twall += tf - t2;
    gaugeData.count += 1;
    gaugeData.matVec += itn;
    if(itn > gaugeData.maxItn) {
        gaugeData.maxItn = itn;
    }
    if(istop == 8) {
        gaugeData.usertol += 1;
    }

    if(istop >= 9) {
        fprintf(stderr, "MINRES returned with istop=%i in %s, exiting.\n",
               istop, __FILE__);
        exit(7);
    }
}

void dynamicClose() {
    if(gaugeData.printDetails > 0){
        printf("dynamicProject:  %i calls (%i usertol), "
               "%i matrix-vector ops, max itn %i\n"
               "                 in %.2f sec (%li wall)\n",
               gaugeData.count, gaugeData.usertol, gaugeData.matVec,
               gaugeData.maxItn,
               gaugeData.tcpu/(float) CLOCKS_PER_SEC, gaugeData.twall);
        printf("dynamicProject:  %.2fvm + %.2fmv = %.2ftotal cpu sec for ops\n",
               gaugeData.tcpu_vm/(float) CLOCKS_PER_SEC,
               gaugeData.tcpu_mv/(float) CLOCKS_PER_SEC,
               (gaugeData.tcpu_vm+gaugeData.tcpu_mv)/(float) CLOCKS_PER_SEC);
        fflush(stdout);
    }

#ifdef USE_MKL
    mkl_free(gaugeData.z);
    mkl_free(gaugeData.b);
    mkl_free(gaugeData.x);
#else
    free(gaugeData.z);
    free(gaugeData.b);
    free(gaugeData.x);
#endif
}

/* Interface for Trlan.
   The extra parameter mvparam is not used in this case. */
void gaugeOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(gaugeData.matrix->rows == abs(nrow));
    assert(mvparam == NULL);
    int k;

    for(k=0; k<ncol; k++) {
        gaugeProduct(&nrow, xin+k*ldx, yout+k*ldy);
    }
}

void gaugeProduct(const integer *vectorLength, const doublereal *x,
                 doublereal *y) {
    clock_t t0 = clock(), t1;
    assert(abs(*vectorLength) == gaugeData.matrix->rows);
    vectorMatrix(gaugeData.matrix, x, gaugeData.z);
    gaugeData.tcpu_vm += (t1=clock()) - t0;
    matrixVector(gaugeData.matrix, gaugeData.z, y);
    gaugeData.tcpu_mv += clock() - t1;
}

/* in and out must be distinct */
void matrixVector(const SparseMatrix *a,
                  const doublereal *in, doublereal *out) {
#ifdef USE_BLOCK
    size_t k;
    double *matp = a->value;
    const integer n = a->blockSize, n2=n*n;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T', normal='N';

    memset(out, 0, a->rows * sizeof(double));
    for(k=0; k<a->blocks; k++, matp+=n2) {
        DGEMV(&trans, &n, &n, &one,
              matp, &n, in + a->j[k], &inc, &one,
              out + a->i[k], &inc);

        if(a->descr == 's' && a->i[k] != a->j[k]) {
            DGEMV(&normal, &n, &n, &one,
                  matp, &n, in + a->i[k], &inc, &one,
                  out + a->j[k], &inc);
        }
    }
#elif defined(USE_MKL)
    sparse_status_t err;
    const double alpha = 1.0, beta=0.0;
    double d;

    if((err=mkl_sparse_d_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
              alpha, a->a, a->descr, in, beta, out, &d)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_dotmv failed %i\n", err);
        exit(69);
    }
#else
    mat_int i, j, k;
    double value;
    memset(out, 0, a->rows * sizeof(double));
    for(k=0; k<a->nonzeros; k++) {
        i = a->i[k];
        j = a->j[k];
        value = a->value[k];
        out[i] += value * in[j];
        if(a->descr == 's' && j != i) {
            out[j] += value * in[i];
        }
    }
#endif
}

void vectorMatrix(const SparseMatrix *a,
                  const doublereal *in, doublereal *out) {
#ifdef USE_BLOCK
    size_t k;
    double *matp = a->value;
    const integer n = a->blockSize, n2 = n*n;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T', normal='N';

    memset(out, 0, a->columns * sizeof(*out));
    for(k=0; k<a->blocks; k++, matp+=n2) {
        DGEMV(&normal, &n, &n, &one,
              matp, &n, in + a->i[k], &inc, &one,
              out + a->j[k], &inc);

        if(a->descr == 's' && a->i[k] != a->j[k]) {
            DGEMV(&trans, &n, &n, &one,
                  matp, &n, in + a->j[k], &inc, &one,
                  out + a->i[k], &inc);
        }
    }
#elif defined(USE_MKL)
    sparse_status_t err;
    const double alpha = 1.0, beta=0.0;
    double d;

    if((err=mkl_sparse_d_dotmv(SPARSE_OPERATION_TRANSPOSE,
              alpha, a->a, a->descr, in, beta, out, &d)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_dotmv failed %i\n", err);
        exit(68);
    }
#else
    mat_int i, j, k;
    double value;
    memset(out, 0, a->columns * sizeof(*out));
    for(k=0; k<a->nonzeros; k++) {
        i = a->i[k];
        j = a->j[k];
        value = a->value[k];
        out[j] += value * in[i];
        if(a->descr == 's' && i != j) {
            out[i] += value * in[j];
        }
    }
#endif
}
