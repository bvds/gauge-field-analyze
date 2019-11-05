/*
    Find shifts in link fields for one step 
    of a saddle-point search.

    Example usage:
    ./shifts ../hess-grad-gauge.json ../shifts.dat
    Valgrind debugging example:
    valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./shifts ../hess-grad-gauge.json junk.dat

      scp hess-grad-gauge.json hess.mtx gauge.mtx grad.dat bvds@192.168.0.35:lattice/gauge-field-analyze/

    Invoke gdb in parallel version:
    mpirun -np 2 xterm -e gdb -ex run --args ./pshifts ../hess-grad-gauge.json junk.out

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> // POSIX
#include <assert.h>
#include <time.h>
#include <omp.h>
#include "shifts.h"
#include "mmio.h"
#ifdef USE_MKL
#include "mkl_spblas.h"
#elif defined(USE_BLIS)
#include "blis/cblas.h"
#elif defined(USE_OPENBLAS)
#include "/opt/OpenBLAS/include/cblas.h"
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


int main(int argc, char **argv){
    char *fdup, *dataPath, *hessFile, *gradFile, *gaugeFile;
    char *options;
    cJSON *jopts, *tmp;
    int i, blockSize, chunkSize, nLargeShifts, wsize, wrank;
    mat_int k, n, partitions, local_n,
        gaugeDimension, local_gaugeDimension;
    FILE *fp;
    int nread, threads;
    long int twall = 0, tcpu = 0, toverall;
    clock_t t1, tt1;
    time_t t2, tt2, tf;

    t1 = clock(); tt1 = t1;
    time(&t2); tt2 = t2;


    /* Initialize MPI */
#ifdef USE_MPI
    if(MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "Failed to initialize MPI.");
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


    /* Read JSON file and use options */
    if(argc <2) {
        fprintf(stderr, "Usage:  ./shifts <parameters.json> <output.dat>\n");
        exit(-1);
    }
    fdup = strdup(argv[1]);
    dataPath = dirname(fdup);
    if(wrank ==0)
        printf("Opening file %s", argv[1]);
    options = readFile(argv[1]);
    jopts = cJSON_Parse(options);
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
    partitions = cJSON_IsNumber(tmp) && tmp->valueint>0?
        (mat_int) tmp->valueint:gaugeDimension/blockSize;
    assert(n%partitions == 0);
    assert((n/partitions)%blockSize == 0);
    assert(gaugeDimension%partitions == 0);
    assert((gaugeDimension/partitions)%blockSize == 0);
    rankSanityTest(partitions);
    local_n = n/partitions*localSize(wrank, wsize, partitions);
    local_gaugeDimension = gaugeDimension/partitions*
        localSize(wrank, wsize, partitions);
    assert(wsize>1 || local_n == n);
    assert(wsize>1 || local_gaugeDimension == gaugeDimension);
    if(wrank ==0)
        printf(" n=%i, gauge=%i blockSize=%i partitions=%i\n",
               n, gaugeDimension, blockSize, partitions);
    // Not used currently
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "chunkSize");
    chunkSize = cJSON_IsNumber(tmp)?tmp->valueint:1;
    assert(chunkSize > 0);

    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "threads");
    threads = cJSON_IsNumber(tmp)?tmp->valueint:1;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(threads);
    if(wrank ==0)
        printf("Setting OpenMP to %i threads\n", threads);
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
    if(wrank ==0)
        printf("Setting MKL to %i threads\n", threads);
    mkl_set_num_threads_local(threads);
#elif defined(USE_BLIS)
    if(wrank ==0)
        printf("Setting Blis to %i threads\n", threads);
    bli_thread_set_num_threads(threads);
#elif defined(USE_OPENBLAS)
    if(wrank ==0)
        printf("Setting OpenBLAS to %i threads\n", threads);
    openblas_set_num_threads(threads);
#else
    assert(threads != 0);  // Suppress compiler warning
#endif


    /*
      Read in Hessian Matrix 
    */
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "hessFile");
    hessFile = catStrings(dataPath, tmp->valuestring);
    SparseMatrix hess;
    sparseMatrixRead(&hess, hessFile, 's', 0, blockSize, partitions,
                     chunkSize, 0, mpicom);
    assert(hess.rows == n);
    assert(hess.columns == n);


    /*
          Read gradient vector
    */
    double *grad = malloc(local_n * sizeof(double)), *gradp = grad, dummy;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "gradFile");
    gradFile = catStrings(dataPath, tmp->valuestring);
    if(wrank ==0)
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
    sparseMatrixRead(&gauge, gaugeFile, 'g', 0, blockSize, partitions,
                     chunkSize, 0, mpicom);
    sparseMatrixRead(&gaugeT, gaugeFile, 'g', 1, blockSize, partitions,
                     chunkSize, 0, mpicom);
    assert(gauge.rows == gaugeDimension);
    assert(gauge.columns == n);
    assert(gaugeT.columns == gaugeDimension);
    assert(gaugeT.rows == n);
    time(&tf);
    tcpu += clock()-t1;
    twall += tf - t2; 

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
                    &nLargeShifts, mpicom);
    linearInit(&hess, local_n, vecs, nLargeShifts, mpicom);
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "linearSolveOptions");
    assert(tmp != NULL);
    linearSolve(local_n, grad, tmp, shifts);
    dynamicClose();


    /* output result */

    t1 = clock();
    time(&t2);
    if(wrank ==0)
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
    time(&tf);
    tcpu += clock()-t1;
    toverall = clock()-tt1;

    sparseMatrixStats(&hess, "hessian");
    sparseMatrixStats(&gauge, "gauge");
    sparseMatrixStats(&gaugeT, "gaugeT");

#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &tcpu, 1, MPI_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &toverall, 1, MPI_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
#endif
    twall += tf - t2;
    if(wrank==0) {
        printf("%s:  input/output in %.2f sec (%li wall)\n",
               __FILE__, tcpu/(float) CLOCKS_PER_SEC, twall);
        printf("%s:  overall time %.2f sec (%li wall)\n",
               __FILE__, (toverall)/(float) CLOCKS_PER_SEC, tf-tt2);
        fflush(stdout);
    }

    sparseMatrixFree(&hess);
    sparseMatrixFree(&gauge);
    sparseMatrixFree(&gaugeT);
    free(fdup);
    free(absVals); free(vecs); free(shifts); free(grad);
    free(hessFile); free(gradFile); free(gaugeFile);
    cJSON_Delete(jopts);
    free(options);
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
