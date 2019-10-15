/*
    Find shifts in link fields for one step 
    of a saddle-point search.

    Example usage:
    ./shifts ../hess-grad-gauge.json ../shifts.dat
    Valgrind debugging example:
    valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./shifts ../hess-grad-gauge.json ../shifts.dat

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
#include "shifts.h"
#include "mmio.h"
#ifdef USE_MKL
#include "mkl_spblas.h"
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
    int i, blockSize, nLargeShifts, wsize, wrank;
    mat_int k, n, local_n, gaugeDimension, local_gaugeDimension;
    FILE *fp;
    MM_typecode matcode;
    int nread;
    mat_int chunkSize;
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
    int threads;
#endif
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
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm mpicom = MPI_COMM_WORLD, *mpicomp = &mpicom;
#else
    void *mpicomp = NULL;
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
    n = cJSON_GetObjectItemCaseSensitive(jopts, "n")->valueint;
    /* Needed for calculation of the cutoff.
       In the parallel case, local size should be
       a multiple of blockSize.  */
    blockSize = cJSON_GetObjectItemCaseSensitive(jopts, "blockSize")->valueint;
    assert(n%blockSize == 0);
    rankSanityTest(n/blockSize);
    local_n = blockSize*localSize(wrank, wsize, n/blockSize);
    assert(wsize>1 || local_n ==n );
    gaugeDimension = cJSON_GetObjectItemCaseSensitive(
                               jopts, "gaugeDimension")->valueint;
    assert(n%gaugeDimension == 0);
    rankSanityTest(gaugeDimension/blockSize);
    local_gaugeDimension = blockSize*localSize(wrank, wsize,
                                               gaugeDimension/blockSize);
    assert(wsize>1 || local_gaugeDimension == gaugeDimension);
    if(wrank ==0)
        printf(" n=%i, blockSize=%i, gauge=%i\n", n, blockSize, gaugeDimension);
    // Not used if using the MKL matrix.
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "chunkSize");
    chunkSize = cJSON_IsNumber(tmp)?tmp->valueint:1;
    assert(chunkSize > 0);

#ifdef USE_MKL
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "threads");
    threads = cJSON_IsNumber(tmp)?tmp->valueint:1;
    if(wrank ==0)
        printf("Setting MKL to %i threads\n", threads);
    mkl_set_num_threads_local(threads);
#endif


    /*
      Read in Hessian Matrix 
    */
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "hessFile");
    hessFile = catStrings(dataPath, tmp->valuestring);
    SparseMatrix hess, *hessp = &hess;
    if(wrank ==0)
        printf("Opening file %s", hessFile);
    if((fp = fopen(hessFile, "r")) == NULL) {
        fprintf(stderr, "%i:  Could not open file %s\n", wrank, hessFile);
        exit(44);
    }
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr, "%i:  Could not process Matrix Market banner.\n",
                wrank);
        exit(701);
    }
    if(!(mm_is_real(matcode) && mm_is_matrix(matcode) && 
         mm_is_coordinate(matcode))) {
        fprintf(stderr, "%i:  Hessian should be real, coordinate.\n", wrank);
        fprintf(stderr, "%i:  Market Market type: [%s]\n", wrank,
                mm_typecode_to_str(matcode));
        exit(702);
    }
    if (mm_read_mtx_crd_size(fp, &(hess.rows), &(hess.columns),
                             &(hess.nonzeros)) !=0) {
        fprintf(stderr, "%i:  Cannot read dimensions\n", wrank);
        exit(703);
    }
    if(wrank ==0)
        printf(" with %i nonzero elements.\n", hess.nonzeros);
    assert(hess.rows == n);
    assert(hess.columns == n);
#ifdef USE_BLOCK
    hess.blockSize = blockSize;
    hess.descr = 's';  // Symmetric matrix, with only one triangle stored.
    sparseMatrixRead(fp, hessp, hessFile, wrank);
    sortMatrix(hessp, (chunkSize+1)/2);
#elif defined(USE_MKL) && !defined(USE_MPI)
    hess.blockSize = blockSize;
    hess.descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    hess.descr.mode = SPARSE_FILL_MODE_LOWER;
    hess.descr.diag = SPARSE_DIAG_NON_UNIT;
    sparseMatrixRead(fp, hessp, hessFile, wrank);
#else
    hess.descr = 's';  // Symmetric matrix, with only one triangle stored.
    sparseMatrixRead(fp, hessp, hessFile, wrank);
    sortMatrix(hessp, (chunkSize+1)/2);
#endif
    fclose(fp);



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
        if(indexRank(k/blockSize, wsize, n/blockSize) == wrank)
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
    SparseMatrix gauge, *gaugep = &gauge;
    if(wrank ==0)
        printf("Opening file %s", gaugeFile);
    fp = fopen(gaugeFile, "r");
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr, "%i:  Could not process Matrix Market banner.\n",
                wrank);
        exit(701);
    }
    if(!(mm_is_real(matcode) && mm_is_matrix(matcode) && 
         mm_is_coordinate(matcode) && mm_is_general(matcode))) {
        fprintf(stderr, "%i:  Gauge matrix should be real, "
                "general, coordinate.\n", wrank);
        fprintf(stderr, "%i:  Market Market type: [%s]\n", wrank,
                mm_typecode_to_str(matcode));
        exit(702);
    }
    if (mm_read_mtx_crd_size(fp, &(gauge.rows), &(gauge.columns),
                             &(gauge.nonzeros)) !=0) {
        fprintf(stderr, "%i:  Cannot read dimensions\n", wrank);
        exit(703);
    }
    if(wrank ==0)
        printf(" with %i nonzero elements.\n", gauge.nonzeros);
    assert(gauge.rows == gaugeDimension);
    assert(gauge.columns == n);
#ifdef USE_BLOCK
    gauge.blockSize = blockSize;
    gauge.descr = 'g';  // General matrix
    sparseMatrixRead(fp, gaugep, gaugeFile, wrank);
    sortMatrix(gaugep, chunkSize);
#elif defined(USE_MKL) && !defined(USE_MPI)
    gauge.blockSize = blockSize;
    gauge.descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    gauge.descr.mode = SPARSE_FILL_MODE_LOWER;
    gauge.descr.diag = SPARSE_DIAG_NON_UNIT;
    sparseMatrixRead(fp, gaugep, gaugeFile, wrank);
#else
    gauge.descr = 'g';  // Symmetric matrix, with only one triangle stored.
    sparseMatrixRead(fp, gaugep, gaugeFile, wrank);
    sortMatrix(gaugep, chunkSize);
#endif
    fclose(fp);
    fflush(stdout);
    time(&tf);
    tcpu += clock()-t1;
    twall += tf - t2; 


    /* Solve it!  */

    double *shifts = malloc(local_n * sizeof(double));
    double *absVals, *vecs;
    int nvals;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "dynamicPartOptions");
    assert(tmp != NULL);
    dynamicInit(local_gaugeDimension, local_n, gaugep, tmp, mpicomp);
    // Debug print:
#if 0
    testOp(hessp, local_n, grad, mpicomp);
#endif
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
    assert(tmp != NULL);
    largeShifts(hessp, tmp, local_n, grad, &absVals, &vecs, &nvals, mpicomp);
    cutoffNullspace(local_n, nvals, jopts, grad, &absVals, &vecs,
                    &nLargeShifts, mpicomp);
    linearInit(hessp, local_n, vecs, nLargeShifts, mpicomp);
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


    sparseMatrixFree(hessp);
    sparseMatrixFree(gaugep);
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
