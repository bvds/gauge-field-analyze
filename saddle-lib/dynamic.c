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
    SparseMatrix *matrixT;
    integer nrow;
    mat_int ncol;
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
    integer itnlim;
    doublereal rtol;
    doublereal *rtolp;
    _MPI_Comm mpicom;
    } gaugeData;


/*
  nrow is the number of rows on this processor
  ncol is the number of columns on this processor
*/
void dynamicInit(const mat_int nrow, const mat_int ncol,
                 SparseMatrix *gauge, SparseMatrix *gaugeT,
                 cJSON *options, _MPI_Comm mpicom) {
    // Let trlan figure out the work array allocation.
    const int lwrk = 0; double *wrk = NULL;
    int i, mev, maxlan, lohi, ned, maxmv, wrank=0;
    double tol = -1.0; // Use default tolerance: sqrt(machine epsilon)
    double *eval, *evec;
    int restart = 3; // restart strategy
    mat_int nonzeros = gauge->blocks*gauge->blockSize*gauge->blockSize;
    trl_info info;
    cJSON *tmp;
#ifdef USE_MPI
    MPI_Comm *mpicomp = &mpicom;
    MPI_Comm_rank(mpicom, &wrank);
#else
    assert(mpicom == NULL);
    void *mpicomp = NULL;
#endif

    gaugeData.matrix = gauge;
    gaugeData.matrixT = gaugeT;
    gaugeData.nrow = nrow;
    gaugeData.ncol = ncol;
    gaugeData.tcpu = 0;
    gaugeData.twall = 0;
    gaugeData.count = 0;
    gaugeData.usertol = 0;
    gaugeData.matVec = 0;
    gaugeData.maxItn = 0;
    gaugeData.printDetails = 0;
    gaugeData.mpicom = mpicom;

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "printDetails");
    if(cJSON_IsNumber(tmp)) {
        gaugeData.printDetails = tmp->valueint;
    } else if(cJSON_IsBool(tmp)) {
        gaugeData.printDetails = cJSON_IsTrue(tmp)?1:0;
    }

    /*
       The gauge configurations are single precision,
       so rtol is normally about 1e-7.
    */
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "rTolerance");
    if(cJSON_IsNumber(tmp)) {
        gaugeData.rtol = tmp->valuedouble;
        gaugeData.rtolp = &gaugeData.rtol;
    } else {
        gaugeData.rtolp = NULL;
    }

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    gaugeData.itnlim = (cJSON_IsNumber(tmp)?tmp->valueint:
                        // Default value appropriate for reorthogonalization
                        (integer) gauge->rows);

    gaugeData.b = MALLOC(nrow*sizeof(*gaugeData.b));
    gaugeData.x = MALLOC(nrow*sizeof(*gaugeData.x));
    gaugeData.z = MALLOC(ncol*sizeof(*gaugeData.z));

    /* Calculate the smallest eigenvalue of gaugeProduct.  
       This will inform the stopping condition for MINRES.

       Tried finding both the smallest and largest eigenvalues,
       but had very slow convergence for large matrices.
    */
    lohi = -1;
    ned = 1;
    // Uses same maxIterations as MINRES.
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxIterations");
    maxmv = cJSON_IsNumber(tmp)?tmp->valueint:(int) nrow*ned;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "maxLanczosVecs");
    maxlan = cJSON_IsNumber(tmp)?tmp->valueint:(int) nrow;

    mev = ned; // Allocate memory for the number of requested eigenpairs
    eval = malloc(mev*sizeof(double));
    evec = malloc(mev*nrow*sizeof(double));
    trl_init_info(&info, nrow, maxlan, lohi, ned, tol, restart, maxmv, mpicomp);
    if(wrank != 0)
        trl_set_debug(&info, 0, "dynamic-trlan-");
    memset(eval, 0, mev*sizeof(double));

    // call TRLAN to compute the eigenvalues
    trlan(gaugeOp, NULL, &info, nrow, mev, eval, evec, nrow, lwrk, wrk);
    if(gaugeData.printDetails > 1) {
        // 2 matrix multiplies, + and * as FLOPS
        trl_print_info(&info, 4*nonzeros);
    } else if(gaugeData.printDetails > 0) {
        trl_terse_info(&info, stdout);
    }
    if(gaugeData.printDetails > 0 && wrank==0) {
        printf("Gauge matrix extreme eigenvalues found:\n");
        for(i=0; i<info.nec; i++) {
            printf("  %le\n", eval[i]);
        }
    }

    if(info.nec>0 && info.stat==0) {
        gaugeData.minEigenvalue = eval[0];
        // gaugeData.maxEigenvalue = eval[info.nec-1];
    } else {
        fprintf(stderr, "%i:  trlan exit with stat=%i, "
                "finding %i of %i eigenpairs.\n",
                wrank, info.stat, info.nec, info.ned);
        exit(8);
    }

    free(eval); free(evec);
}

void dynamicProject(const integer n, double *v, double *normDiff) {
    doublereal shift = 0.0, *maxxnormp = NULL;
    // Explicit value for these, so we can force MINRES alogrithm.
    doublereal trancond, acondlim = 1.0e15;
    logical *disablep = NULL;
    integer nout, *noutp = NULL, istop, itn;
    doublereal rnorm, arnorm, xnorm, anorm, acond,
        abstol, *abstolp = NULL;
    doublereal normv;
    const integer one=1;
    const doublereal minusone = -1.0;
    int wrank;
#ifdef USE_MPI
    MPI_Comm mpicom = gaugeData.mpicom;
    MPI_Comm_rank(mpicom, &wrank);
#else
    wrank = 0;
#endif
    clock_t t1;
    time_t t2, tf;

    assert(gaugeData.ncol == (mat_int) n);

    t1 = clock();
    time(&t2);

    /*
      Tolerance of the projection relative to the initial vector.
      For now, tie this to rTolerance.  May want separate flag?
    */
    if(gaugeData.rtolp != NULL) {
        normv = DNRM2(&n, v, &one);
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &normv, 1, MPI_DOUBLE_PRECISION,
                      MPI_SUM, mpicom);
#endif
        abstol = gaugeData.rtol*normv*sqrt(gaugeData.minEigenvalue);
        abstolp = &abstol;
#if 0
        if(wrank == 0)
            printf("Setting abstolp=%le normin=%le val=%le sqrt=%le\n",
                   abstol, normin, gaugeData.minEigenvalue,
                   sqrt(gaugeData.minEigenvalue));
#endif
    }
    if(gaugeData.printDetails > 2) {
        nout = 6;  // Write to stdout
        noutp = &nout;
    }

    matrixVector(gaugeData.matrix, n, v,
                 gaugeData.nrow, gaugeData.b);

    trancond = acondlim; // Always use MINRES
    MINRESQLP(&gaugeData.nrow, gaugeProduct, gaugeData.b,
              &shift, NULL, NULL, disablep, noutp, &gaugeData.itnlim,
              gaugeData.rtolp, abstolp, maxxnormp, &trancond, &acondlim,
              gaugeData.x, &istop, &itn, &rnorm, &arnorm,
              &xnorm, &anorm, &acond);

    matrixVector(gaugeData.matrixT,
                 gaugeData.nrow, gaugeData.x,
                 n, gaugeData.z);
    // This could be combined with the above matrix product, BLAS style.
    DAXPY(&n, &minusone, gaugeData.z, &one, v, &one);
    if(normDiff != NULL){
        *normDiff = DNRM2(&n, gaugeData.z, &one);
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, normDiff, 1, MPI_DOUBLE,
                      MPI_SUM, mpicom);
#endif
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
        fprintf(stderr, "%i: MINRES returned with istop=%i in %s, exiting.\n",
                wrank, istop, __FILE__);
        exit(7);
    }
}

void dynamicClose() {
    int wrank;
#ifdef USE_MPI
    MPI_Comm mpicom = gaugeData.mpicom;
    MPI_Comm_rank(mpicom, &wrank);
    MPI_Allreduce(MPI_IN_PLACE, &gaugeData.tcpu, 1, MPI_LONG,
                  MPI_SUM, mpicom);
    MPI_Allreduce(MPI_IN_PLACE, &gaugeData.tcpu_vm, 1, MPI_LONG,
                  MPI_SUM, mpicom);
    MPI_Allreduce(MPI_IN_PLACE, &gaugeData.tcpu_mv, 1, MPI_LONG,
                  MPI_SUM, mpicom);
#else
    wrank = 0;
#endif
    if(gaugeData.printDetails > 0 && wrank == 0){
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

    FREE(gaugeData.z);
    FREE(gaugeData.b);
    FREE(gaugeData.x);
}

/* Interface for Trlan.
   The extra parameter mvparam is not used in this case. */
void gaugeOp(const int nrow, const int ncol, const double *xin, const int ldx,
             double *yout, const int ldy, void* mvparam) {
    assert(mvparam == NULL);
    int k;
    integer n = nrow;

    for(k=0; k<ncol; k++) {
        gaugeProduct(&n, xin+k*ldx, yout+k*ldy);
    }
}

void gaugeProduct(const integer *vectorLength, const doublereal *x,
                  doublereal *y) {
    clock_t t0 = clock(), t1;
    mat_int vl = *vectorLength;

    assert(gaugeData.nrow == *vectorLength);

    matrixVector(gaugeData.matrixT, vl, x,
                 gaugeData.ncol, gaugeData.z);
    gaugeData.tcpu_vm += (t1=clock()) - t0;
    matrixVector(gaugeData.matrix,
                 gaugeData.ncol, gaugeData.z, vl, y);
    gaugeData.tcpu_mv += clock() - t1;
}
