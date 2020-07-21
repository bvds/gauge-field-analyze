#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "shifts.h"
/*
     Apply cutoff to eigenpairs.
     The norm is an L-2 norm applied to each nc^2-1 rows
     with an L-infinity norm applied to that.

     This differs from the norm defined in the paper by
     a factor of sqrt(2), but eigenCutoffMax, eigenCutoff2, 
     and eigenCutoffRescale should match the associated lambda
     parameters in the paper.

     Returns filtered vals & vecs arrays.
     The associated memory is re-allocated.

     n is the vector length for the local process.
*/
void cutoffNullspace(mat_int n, int nvals, cJSON *options,
                     double *grad,
                     double **vals, double **vecs, int *nLargeShifts,
                     cJSON *jout, _MPI_Comm mpicom) {
    mat_int j, blocks;
    int i, wrank;
    int testConcave, testMax, test2,
        countConcave = 0, countMax = 0, count2 = 0;
    int na;
    int firstValue = -1;
    int lastValue = -1;
    double norm2, vecnorm2,  maxnorm2, vecdotgrad;
    double cutoffMax, cutoff2, zzz;
    cJSON *tmp;
    cJSON_bool removeConcave, squareHess;
#ifdef USE_MPI
    MPI_Comm_rank(mpicom, &wrank);
#else
    assert(mpicom == NULL);
    wrank = 0;
#endif

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "blockSize");
    assert(cJSON_IsNumber(tmp));
    na = tmp->valueint;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "removeConcave");
    removeConcave = cJSON_IsBool(tmp)?cJSON_IsTrue(tmp):1==0;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenCutoffMax");
    cutoffMax = cJSON_IsNumber(tmp)?tmp->valuedouble:1.0;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenCutoff2");
    cutoff2 = cJSON_IsNumber(tmp)?tmp->valuedouble:1.0;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "eigenCutoffRescale");
    zzz = cJSON_IsNumber(tmp)?tmp->valuedouble:1.0;

    // Sanity test for squaring the Hessian
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "largeShiftOptions");
    assert(tmp != NULL);
    tmp  = cJSON_GetObjectItemCaseSensitive(tmp, "squareHess");
    squareHess = cJSON_IsBool(tmp)?cJSON_IsTrue(tmp):1==1;
    assert(removeConcave != squareHess);

    // Color blocks should be evenly divided among processes
    assert(n%na == 0);

    *nLargeShifts = 0;
#if 0  // Debug print
    printf("cutoffNullSpace for n=%i, nvals=%i\nEigenvalues:\n", n, nvals);
    for(i=0; i<nvals; i++) {
        printf("  %le\n", vals[i]);
    }
#endif
    for(i=0; i<nvals; i++) {
        maxnorm2 = 0.0;
        vecnorm2 = 0.0;
        norm2 = 0.0;
        vecdotgrad = 0.0;
        for(j=0; j<n; j++) {
            norm2 += pow((*vecs)[i*n + j], 2);
            vecdotgrad += (*vecs)[i*n + j] * grad[j];
            if((j+1)%na == 0) {
                if(norm2 > maxnorm2) {
                    maxnorm2 = norm2;
                }
                vecnorm2 += norm2;
                norm2 = 0;
            }
        }
        blocks = n/na;
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &vecnorm2, 1, MPI_DOUBLE,
                      MPI_SUM, mpicom);
        MPI_Allreduce(MPI_IN_PLACE, &vecdotgrad, 1, MPI_DOUBLE,
                      MPI_SUM, mpicom);
        MPI_Allreduce(MPI_IN_PLACE, &maxnorm2, 1, MPI_DOUBLE,
                      MPI_MAX, mpicom);
        MPI_Allreduce(MPI_IN_PLACE, &blocks, 1, _MPI_MAT_INT,
                      MPI_SUM, mpicom);
#endif
        testConcave = removeConcave && (*vals)[i] < 0.0;
        if(testConcave)
            countConcave++;
        /* If eigenvectors are normalized, then vecnorm2
           is not needed. */
        testMax = cutoffMax*zzz*fabs((*vals)[i])*vecnorm2 <=
            fabs(vecdotgrad)*sqrt(maxnorm2);
        if(testMax)
            countMax++;
        test2 = cutoff2*zzz*fabs((*vals)[i])*sqrt(vecnorm2*blocks) <=
            fabs(vecdotgrad);
        if(test2)
            count2++;
        if(testMax || test2 || testConcave) {
            if(firstValue < 0)
                firstValue = i;
            lastValue = i;
            if(*nLargeShifts < i) {
                (*vals)[*nLargeShifts] = (*vals)[i];
                memcpy(*vecs+(*nLargeShifts)*n, *vecs+i*n, n*sizeof(double));
            }
            *nLargeShifts += 1;
        }
    }
    *vals = realloc(*vals, (*nLargeShifts)*sizeof(double));
    *vecs = realloc(*vecs, n*(*nLargeShifts)*sizeof(double));

    cJSON_AddNumberToObject(jout, "firstValue", firstValue);
    cJSON_AddNumberToObject(jout, "lastValue", lastValue);
    cJSON_AddNumberToObject(jout, "nvals", nvals);
    cJSON_AddNumberToObject(jout, "nvalsMax", countMax);
    cJSON_AddNumberToObject(jout, "nvals2", count2);
    cJSON_AddNumberToObject(jout, "nvalsConcave", countConcave);
    if(wrank == 0) {
        printf("cutoffNullSpace:  %u (max %i, norm %i, concave %i) "
               "of %u eigenpairs between [%i, %i]\n",
               *nLargeShifts, countMax, count2, countConcave,
               nvals, firstValue, lastValue);
    }
}
