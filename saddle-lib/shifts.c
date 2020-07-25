/*
    Find shifts in link fields for one step 
    of a saddle-point search.

    Example usage:
    ./shifts hess-grad-gauge.json shifts.dat out.json
    Valgrind debugging example:
    valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./shifts ../hess-grad-gauge.json junk.dat junk.json

    Invoke gdb in parallel version:
    mpirun -np 2 xterm -e gdb -ex run --args ./pshifts ../hess-grad-gauge.json junk.out junk.json

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> // POSIX
#include <assert.h>
#include <omp.h>
#include "shifts.h"
#include "mmio.h"
#ifdef USE_MKL
#include "mkl_spblas.h"
#elif defined(USE_BLIS)
#include "blis/cblas.h"
#elif defined(USE_OPENBLAS)
#include "cblas.h"
#endif


char *catStrings(const char *str1, const char *str2);

char *catStrings(const char *str1, const char *str2) {
    const char d[] = "/";
    // Add a byte for the string terminator.
    char *x = malloc((strlen(str1) + strlen(d) + strlen(str2) + 1)*
                     sizeof(char));
    char *y = stpcpy(x, str1);
    y = stpcpy(y, d);
    stpcpy(y, str2);
    // printf("catStrings %s+%s -> %s\n", str1, str2, x);
    return x;
}

char *readFile(char *filename) {
    int k;
    FILE *f = fopen(filename, "rt");
    assert(f != NULL);
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buffer = (char *) malloc(length + 1);
    buffer[length] = '\0';
    k = fread(buffer, 1, length, f);
    assert(k == length);
    fclose(f);
    return buffer;
}

int getPrintDetails(cJSON *jopts, int defaultLevel) {
    cJSON *tmp = cJSON_GetObjectItemCaseSensitive(jopts, "printDetails");
    if(cJSON_IsNumber(tmp))
        return tmp->valueint;
    else if(cJSON_IsBool(tmp))
        return cJSON_IsTrue(tmp)?defaultLevel:0;
    else
        return defaultLevel;
}

#ifndef USE_MPI
double tdiff(struct timeval t1, struct timeval t0) {
    return 1.0e-6 * (double) (t1.tv_usec - t0.tv_usec) +
        (double) (t1.tv_sec - t0.tv_sec);
}
#endif

int main(int argc, char **argv){
    char *fdup, *dataPath, *hessFile, *gradFile, *gaugeFile;
    char *options;
    cJSON *jopts, *tmp, *jout;
    int i, blockSize, chunkSize, nLargeShifts, wsize, wrank;
    mat_int k, n, partitions, local_n,
        gaugePartitions, gaugeDimension, local_gaugeDimension;
    FILE *fp;
    int nread, threads, printDetails;
    TIME_TYPE t0, t1, t2;
    double tio = 0.0, tcalc = 0.0;

    jout = cJSON_CreateObject();

    /* Initialize MPI */
#ifdef USE_MPI
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &i);
    if(i != MPI_THREAD_FUNNELED) {
        fprintf(stderr, "Warning MPI did not provide MPI_THREAD_FUNNELED\n");
        exit(321);
    }
    MPI_Comm mpicom = MPI_COMM_WORLD;
    MPI_Comm_size(mpicom, &wsize);
    MPI_Comm_rank(mpicom, &wrank);
#else
    void *mpicom = NULL;
    wsize = 1;
    wrank = 0;
#endif

    SET_TIME(t0);

    /* Read JSON file and use options */
    if(argc <2) {
        fprintf(stderr, "Usage:  ./shifts <parameters.json> <output.dat>\n");
        exit(-1);
    }
    fdup = strdup(argv[1]);
    dataPath = dirname(fdup);
    options = readFile(argv[1]);
    jopts = cJSON_Parse(options);
    printDetails = getPrintDetails(jopts, 3);
    /* Needed for calculation of the cutoff.
       In the parallel case, local size should be
       a multiple of blockSize.  */
    blockSize = cJSON_GetObjectItemCaseSensitive(jopts, "blockSize")->valueint;
    n = cJSON_GetObjectItemCaseSensitive(jopts, "n")->valueint;
    assert(n%blockSize == 0);
    gaugeDimension = cJSON_GetObjectItemCaseSensitive(
                               jopts, "gaugeDimension")->valueint;
    assert(gaugeDimension%blockSize == 0);
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "partitions");
    if(cJSON_IsNumber(tmp)) {
        partitions = (mat_int) tmp->valueint;
        gaugePartitions = partitions;
    } else {
        partitions = n/blockSize;
        gaugePartitions = gaugeDimension/blockSize;
    }
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "gaugePartitions");
    if(cJSON_IsNumber(tmp)) {
        gaugePartitions = (mat_int) tmp->valueint;
    }
    assert(n%partitions == 0);
    assert((n/partitions)%blockSize == 0);
    assert(gaugeDimension%gaugePartitions == 0);
    assert((gaugeDimension/gaugePartitions)%blockSize == 0);
    rankSanityTest(partitions);
    rankSanityTest(gaugePartitions);
    local_n = n/partitions*localSize(wrank, wsize, partitions);
    local_gaugeDimension = gaugeDimension/gaugePartitions*
        localSize(wrank, wsize, gaugePartitions);
    assert(wsize>1 || local_n == n);
    assert(wsize>1 || local_gaugeDimension == gaugeDimension);
    if(wrank == 0 && printDetails > 2)
        printf("Parameter file %s:  n=%i gauge=%i "
               "blockSize=%i partitions=%i gaugePartitions = %i\n",
               argv[1], n, gaugeDimension, blockSize, partitions,
               gaugePartitions);
    // Not used currently
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "chunkSize");
    chunkSize = cJSON_IsNumber(tmp)?tmp->valueint:1;
    assert(chunkSize > 0);

    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "threads");
    threads = cJSON_IsNumber(tmp)?tmp->valueint:1;
#ifdef USE_MKL
    /*
      On the T480 laptop, calculations are limited by memory
      latency/bandwidth:  too many threads causes thread 
      contention.

      8^3 lattice, nc=3 total execution time:
      Threads    cpu time   wall time
      1          884s       884s
      2          1362s      682s
      3          2042s      683s
      4          2878s      722s
    */
    mkl_set_num_threads_local(threads);
    char libName[] = " & MKL";
#elif defined(USE_BLIS)
    bli_thread_set_num_threads(threads);
    char libName[] = " & Blis";
#elif defined(USE_OPENBLAS)
    openblas_set_num_threads(threads);
    char libName[] = " & OpenBLAS";
#else
    char libName[] = "";
#endif
    /* Tried disabling dynamic teams,
       but this degraded performance. */
    omp_set_num_threads(threads);
    if(wrank == 0 && printDetails > 2)
        printf("Setting OpenMP%s to %i threads\n",
               libName, threads);


    /*
      Read in Hessian Matrix 
    */
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "hessFile");
    hessFile = catStrings(dataPath, tmp->valuestring);
    SparseMatrix hess;
    sparseMatrixRead(&hess, hessFile, 's', 0, blockSize, partitions,
                     chunkSize, printDetails, mpicom);
    assert(hess.rows == n);
    assert(hess.columns == n);


    /*
          Read gradient vector
    */
    double *grad = malloc(local_n * sizeof(double)), *gradp = grad, dummy;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "gradFile");
    gradFile = catStrings(dataPath, tmp->valuestring);
    if(wrank == 0 && printDetails > 2)
        printf("Opening file %s for a vector of length %i\n", gradFile, n);
    fp = fopen(gradFile, "r");
    for(k=0; k<n; k++){
        if(indexRank(k/(n/partitions), wsize, partitions) == wrank)
            nread = fscanf(fp, "%le", gradp++);
        else
            nread = fscanf(fp, "%le", &dummy);
        if(nread < 1) {
            fprintf(stderr, "%i:  Error reading %s on line %i\n",
                    wrank, gradFile, k);
            break;
        }
    }
    assert(local_n == gradp - grad);
    fclose(fp);


    /*
           Read gauge-invariant shifts matrix
    */
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "gaugeFile");
    gaugeFile = catStrings(dataPath, tmp->valuestring);
    SparseMatrix gauge, gaugeT;
    sparseMatrixRead(&gauge, gaugeFile, 'g', 0, blockSize, gaugePartitions,
                     chunkSize, printDetails, mpicom);
    sparseMatrixRead(&gaugeT, gaugeFile, 'g', 1, blockSize, partitions,
                     chunkSize, printDetails, mpicom);
    assert(gauge.rows == gaugeDimension);
    assert(gauge.columns == n);
    assert(gaugeT.columns == gaugeDimension);
    assert(gaugeT.rows == n);

    ADD_TIME(tio, t1, t0);
#if 0 // Debug:  Test matrix-vector multiplication, printing result  
    testMatrixVector(&gauge, grad);
    exit(0);
#endif

    /* Solve it!  */

    double *shifts = malloc(local_n * sizeof(double));
    double *absVals, *vecs;
    int nvals;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "dynamicPartOptions");
    assert(tmp != NULL);
    dynamicInit(local_gaugeDimension, local_n, &gauge, &gaugeT, tmp, mpicom);
    // Debug:  Test hessOp, printing result.
#if 0
    testOp(&hess, local_n, grad, mpicom);
#endif
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
    assert(tmp != NULL);
    largeShifts(&hess, tmp, local_n, grad, &absVals, &vecs, &nvals, mpicom);
    cutoffNullspace(local_n, nvals, jopts, grad, &absVals, &vecs,
                    &nLargeShifts, jout, mpicom);
    linearInit(&hess, local_n, vecs, nLargeShifts, mpicom);
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "linearSolveOptions");
    assert(tmp != NULL);
    linearSolve(local_n, grad, tmp, shifts);
    dynamicClose();


    /* output result */

    ADD_TIME(tcalc, t2, t1);
    if(wrank == 0 && printDetails>2)
        printf("Opening output file %s\n", argv[2]);
    for(i=0; i<wsize; i++) {
        if(wrank == i) {
            fp = fopen(argv[2], wrank==0?"w":"a");
            for(k=0; k<local_n; k++) {
                fprintf(fp, "%.17e\n", shifts[k]);
            }
            fclose(fp);
        }
#ifdef USE_MPI
        MPI_Barrier(mpicom);
#endif
    }
    if(wrank == 0) {
        char *joutString = cJSON_Print(jout);
        if(printDetails>2)
            printf("Opening output file %s\n", argv[2]);
        fp = fopen(argv[3], "w");
        fprintf(fp, "%s\n", joutString);
        fclose(fp);
        free(joutString);
    }

    if(printDetails > 1) {
        sparseMatrixStats(&hess, "hessian");
        sparseMatrixStats(&gauge, "gauge");
        sparseMatrixStats(&gaugeT, "gaugeT");
    }

    ADD_TIME(tio, t0, t2);
#ifdef USE_MPI
    /* This averaging isn't really necessary ...
       It might make more sense to find the maximum. */
    MPI_Allreduce(MPI_IN_PLACE, &tio, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &tcalc, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    tio /= wsize;
    tcalc /= wsize;
#endif
    if(wrank==0 && printDetails > 0)
        printf("%s:  init and i/o %.2f s, calculation %.2f s\n",
               __FILE__, tio, tcalc);

    sparseMatrixFree(&hess);
    sparseMatrixFree(&gauge);
    sparseMatrixFree(&gaugeT);
    free(fdup);
    free(absVals); free(vecs); free(shifts); free(grad);
    free(hessFile); free(gradFile); free(gaugeFile);
    cJSON_Delete(jopts);
    free(options);
    cJSON_Delete(jout);
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
