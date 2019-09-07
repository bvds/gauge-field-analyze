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
  the shifts using MINRES/MINRES-QLP
*/

struct {
    SparseMatrix *matrix;
    cJSON *options;
    doublereal *z;
    doublereal *x;
    doublereal *b;
    unsigned long int tcpu;
    unsigned long int twall;
    unsigned int count;
    unsigned int usertol;
    unsigned int matVec;
    int printDetails;
    doublereal minEigenvalue;
    doublereal maxEigenvalue;
} gaugeData;



void dynamicInit(SparseMatrix *gauge, cJSON *options) {
    // Let trlan figure out the work array allocation.
    const int lwrk = 0; double *wrk = NULL;
    int mev, maxlan, lohi, ned, maxmv;
    double tol = -1.0; // Use default tolerance: sqrt(machine epsilon)
    double *eval, *evec;
    int nrow; // number of rows on this processor
    int restart = 7; // restart strategy
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
    gaugeData.printDetails = 0;

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "printDetails");
    if(cJSON_IsNumber(tmp)) {
        gaugeData.printDetails = tmp->valueint;
    } else if(cJSON_IsBool(tmp)) {
        gaugeData.printDetails = cJSON_IsTrue(tmp)?1:0;
    }

    gaugeData.b = malloc(rows(gauge)*sizeof(double));
    gaugeData.x = malloc(rows(gauge)*sizeof(double));
    gaugeData.z = malloc(columns(gauge)*sizeof(doublereal));


    /* Calculate the largest and smallest eigenvalues
       of gaugeProduct.  These will inform the stopping condition
       for MINRES. */
    nrow = rows(gauge);
    ned = 4;  // In principle, 2 should be enough.
    lohi = 0;
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "maxIterations");
    maxmv = cJSON_IsNumber(tmp)?tmp->valueint:ned*nrow;
    tmp  = cJSON_GetObjectItemCaseSensitive(gaugeData.options, "maxLanczosVecs");
    maxlan = cJSON_IsNumber(tmp)?tmp->valueint:nrow;

    mev = ned; // Allocate memory for the number of requested eigenpairs
    eval = malloc(mev*sizeof(double));
    evec = malloc(mev*nrow*sizeof(double));
    trl_init_info(&info, nrow, maxlan, lohi, ned, tol, restart, maxmv, 0);
    memset(eval, 0, mev*sizeof(double));

    // call TRLAN to compute the eigenvalues
    trlan(gaugeOp, NULL, &info, nrow, mev, eval, evec, nrow, lwrk, wrk);
    if(gaugeData.printDetails > 1) {
        // 2 matrix multplies, with 1 add and 1 multiply per nonzero element.
        trl_print_info(&info, 4*nonzeros(gauge));
    } else if(gaugeData.printDetails > 0) {
        trl_terse_info(&info, stdout);
    }
    if(gaugeData.printDetails > 0) {
        printf("Gauge matrix extreme eigenvalues found:\n");
        for(i=0; i<info.nec; i++) {
            printf("  %le\n", eval[i]);
        }
    }

    if(info.nec>=2 && info.stat==0) {
        gaugeData.minEigenvalue = eval[0];
        gaugeData.maxEigenvalue = eval[info.nec-1];
    } else {
        printf("trlan exit with stat=%i, finding %i of %i eigenpairs.\n",
               info.stat, info.nec, info.ned);
        exit(8);
    }

    free(eval); free(evec);
}

/* "in" and "out" may overlap */
void dynamicProject(const int n, double *in, double *out) {
    int i, rr;
    doublereal shift = 0.0, *maxxnormp = NULL;
    // Explicit value for these, so we can force MINRES alogrithm.
    doublereal trancond, acondlim = 1.0e15;
    logical *disablep = NULL;
    integer nout, *noutp = NULL, istop, itn, itnlim;
    cJSON *tmp;
    doublereal rnorm, arnorm, xnorm, anorm, acond,
        rtol, *rtolp = NULL, abstol, *abstolp = NULL;
    doublereal normin;
    const integer one=1;
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
        normin = DNRM2(&n, in, &one);
        abstol = rtol*normin*sqrt(gaugeData.minEigenvalue);
        abstolp = &abstol;
#if 0
        printf("Setting abstolp=%le normin=%le val=%le sqrt=%le\n",
               abstol, normin, gaugeData.minEigenvalue, sqrt(gaugeData.minEigenvalue));
#endif
    }
    if(gaugeData.printDetails > 2) {
        nout = 6;  // Write to stdout
        noutp = &nout;
    }

    // Sanity test for gaugeData.z
    assert(n == columns(gaugeData.matrix));

    matrixVector(gaugeData.matrix, in, gaugeData.b);
  
    trancond = acondlim; // Always use MINRES
    rr = rows(gaugeData.matrix);
    MINRESQLP(&rr, gaugeProduct, gaugeData.b,
              &shift, NULL, NULL, disablep, noutp, &itnlim,
              rtolp, abstolp, maxxnormp, &trancond, &acondlim,
              gaugeData.x, &istop, &itn, &rnorm, &arnorm,
              &xnorm, &anorm, &acond);

    vectorMatrix(gaugeData.matrix, gaugeData.x, gaugeData.z);
    // This could be combined with the matrix product, BLAS style.
    for(i=0; i<n; i++) {
        out[i] = in[i] - gaugeData.z[i];
    }

    time(&tf);
    gaugeData.tcpu += clock()-t1;
    gaugeData.twall += tf - t2; 
    gaugeData.count += 1;
    gaugeData.matVec += itn;
    if(istop == 8) {
        gaugeData.usertol += 1;
    }

    if(istop >= 9) {
        printf("MINRES returned with istop=%i in %s, exiting.\n",
               istop, __FILE__);
        exit(7);
    }
}

void dynamicClose() {
    if(gaugeData.printDetails > 0){
        printf("dynamicProject:  %i calls (%i usertol), "
               "%i matrix-vector ops,\n"
               "                 in %.2f sec (%li wall)\n",
               gaugeData.count, gaugeData.usertol, gaugeData.matVec, 
               gaugeData.tcpu/(float) CLOCKS_PER_SEC, gaugeData.twall);
    }

    free(gaugeData.z); gaugeData.z = NULL;
    free(gaugeData.b);
    free(gaugeData.x);
}

/* Interface for Trlan.
   The extra parameter mvparam is not used in this case. */
void gaugeOp(const int nrow, const int ncol, const double *xin, const int ldx,
	    double *yout, const int ldy, void* mvparam) {
    assert(rows(gaugeData.matrix) == nrow);
    assert(mvparam == NULL);
    int k;

    for(k=0; k<ncol; k++) {
        gaugeProduct(&nrow, xin+k*ldx, yout+k*ldy);
    }
}

int gaugeProduct(const integer *vectorLength, const doublereal *x,
                 doublereal *y) {
    assert(*vectorLength == rows(gaugeData.matrix));
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
