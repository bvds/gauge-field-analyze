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

     Returns filtered vals & vecs arrays.
     The associated memory is re-allocated.

     n is the vector length for the local process.
*/
void cutoffNullspace(mat_int n, int nvals, cJSON *options,
                     double *grad,
                     double **vals, double **vecs, int *nLargeShifts,
                     cJSON *jout, _MPI_Comm mpicom) {
    mat_int j;
    int i, wrank;
    int na;
    int firstValue = -1;
    int lastValue = -1;
    double norm2, vecnorm2,  maxnorm2, vecdotgrad;
    double cutoff, zzz;
    cJSON *tmp;
#ifdef USE_MPI
    MPI_Comm_rank(mpicom, &wrank);
#else
    assert(mpicom == NULL);
    wrank = 0;
#endif

    tmp  = cJSON_GetObjectItemCaseSensitive(options, "blockSize");
    assert(cJSON_IsNumber(tmp));
    na = tmp->valueint;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "largeShiftCutoff");
    cutoff = cJSON_IsNumber(tmp)?tmp->valuedouble:1.0;
    tmp  = cJSON_GetObjectItemCaseSensitive(options, "removeCutoff");
    zzz = cJSON_IsNumber(tmp)?tmp->valuedouble:1.0;

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
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &vecnorm2, 1, MPI_DOUBLE,
                      MPI_SUM, mpicom);
        MPI_Allreduce(MPI_IN_PLACE, &vecdotgrad, 1, MPI_DOUBLE,
                      MPI_SUM, mpicom);
        MPI_Allreduce(MPI_IN_PLACE, &maxnorm2, 1, MPI_DOUBLE,
                      MPI_MAX, mpicom);
#endif
        /* If eigenvectors are normalized, then vecnorm2
           is not needed. */
        if(cutoff * zzz * fabs((*vals)[i]) * vecnorm2 <=
           fabs(vecdotgrad) * sqrt(maxnorm2)) {
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
    if(wrank == 0) {
        printf("cutoffNullSpace:  %u of %u zeros between [%i, %i]\n",
               *nLargeShifts, nvals, firstValue, lastValue);
    }
}
