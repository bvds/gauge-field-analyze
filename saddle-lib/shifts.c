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
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "shifts.h"
#include "mmio.h"
#ifdef USE_MKL
#include "mkl_spblas.h"
#endif

/*
  Divide vectors and matrices up for mpi.  
  Without mpi, wrank = 0 and wsize = 1 gives
  the desired single-process behavior.
 */
int indexRank(const mat_int i, const int wsize, const mat_int n) {
    if(i/(n/wsize + 1) < n%wsize) 
        return i/(n/wsize + 1);
    else
        return (i-n%wsize)/(n/wsize);
}
mat_int localSize(const unsigned int wrank, const int wsize, const mat_int n) {
    return n/wsize + ((wrank<(n%wsize))?1:0);
}
mat_int rankIndex(const unsigned int wrank, const int wsize, const mat_int n) {
    if(wrank < n%wsize)
        return wrank*(n/wsize + 1);
    else
        return wrank*(n/wsize) + n%wsize;
}

// Sanity tests for localSize, rankIndex, and indexRank
void rankSanityTest(mat_int n) {
    int i, wsize = 1;
    mat_int k;
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
#endif
    k = 0;
    for(i=0; i<wsize; i++) {
        assert(i == indexRank(rankIndex(i, wsize, n), wsize, n));
        assert(localSize(i, wsize, n) == rankIndex(i + 1, wsize, n)
               - rankIndex(i, wsize, n));
        k += localSize(i, wsize, n);
    }
    assert(k == n);
}

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

#ifdef USE_BLOCK
void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName, int wrank);
void blockFree(SparseMatrix *mat);

void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName, int wrank) {
    mat_int j;
    mat_int i, k, ii, jj, lastii = 0, blockCols = 0;
    mat_int *blockColp;
    const mat_int block = mat->blockSize,
        maxBlockCol = mat->columns/block;
    double value, *blockRowValue, *blockValuep;
    int *blockColFlag, nread;
    assert(mat->rows%block == 0);
    assert(mat->columns%block == 0);
    mat->blocks = 0;
    blockRowValue = malloc(mat->columns*block*sizeof(*blockRowValue));
    blockColFlag = malloc(maxBlockCol*sizeof(*blockColFlag));
    blockColp = malloc(maxBlockCol*sizeof(*blockColp));
#ifdef USE_MKL
    mat->value = mkl_malloc(0, MALLOC_ALIGN);
    mat->i = mkl_malloc(0, MALLOC_ALIGN);
    mat->j = mkl_malloc(0, MALLOC_ALIGN);
#else
    mat->value = NULL;
    mat->i = NULL;
    mat->j = NULL;
#endif
    memset(blockRowValue, 0, mat->columns*block*sizeof(double));
    memset(blockColFlag, 0, maxBlockCol*sizeof(*blockColFlag));
    for(j=0;; j++){
        if(j<mat->nonzeros) {
            nread = fscanf(fp, "%u%u%le", &ii, &jj, &value);
            if(nread < 3) {
                fprintf(stderr, "%i:  Error reading %s, element %i\n",
                        wrank, fileName, j);
                break;
            }
            ii -= 1; jj -= 1; // switch to zero-based indexing
            assert(ii>=lastii);
        }
        if(j==mat->nonzeros || (lastii<ii && ii%block == 0)) {
            // Retire block row and reset for a new one
#ifdef USE_MKL
            mat->value = mkl_realloc(mat->value,
                     (mat->blocks + blockCols)*block*block*sizeof(double));
            mat->i = mkl_realloc(mat->i,
                     (mat->blocks + blockCols)*sizeof(*mat->i));
            mat->j = mkl_realloc(mat->j,
                     (mat->blocks + blockCols)*sizeof(*mat->j));
#else
            mat->value = realloc(mat->value,
                    (mat->blocks + blockCols)*block*block*sizeof(double));
            mat->i = realloc(mat->i,
                    (mat->blocks + blockCols)*sizeof(*mat->i));
            mat->j = realloc(mat->j,
                    (mat->blocks + blockCols)*sizeof(*mat->j));
#endif
            for(i=0; i<blockCols; i++) {
                k = blockColp[i];
#if 0
                printf("  %i %i mat->blocks=%i blockCols=%i\n",
                       j, k, mat->blocks, blockCols);
#endif
                blockValuep = blockRowValue+k*block*block;
                memcpy(mat->value+block*block*mat->blocks,
                       blockValuep, block*block*sizeof(double));
                memset(blockValuep, 0, block*block*sizeof(double));
                mat->i[mat->blocks] = (lastii/block)*block;
                mat->j[mat->blocks] = k*block;
                mat->blocks += 1;
                blockColFlag[k] = 0;
            }
            blockCols = 0;
        }
        if(j==mat->nonzeros)
            break;
        lastii = ii;
        // Mark this block
        k = jj/block;
        if(!blockColFlag[k]) {
            blockColp[blockCols] = k;
            blockCols += 1;
            blockColFlag[k] = 1;
        }
        // Add value to this block
        blockRowValue[(k*block + ii%block)*block + jj%block] = value;
        /* fill diagonal blocks of symmetric matrix */
        if(mat->descr == 's' && ii/block == k && ii!=jj) {
            blockRowValue[(k*block + jj%block)*block + ii%block] = value;
        }
    }
    free(blockRowValue);
    free(blockColFlag);
    free(blockColp);
#if 0
    printf("Print matrix for file %s\n", fileName);
    for(k=0; k<mat->blocks; k++) {
        printf("%i %i  ", mat->i[k] + 1, mat->j[k] + 1);
        for(ii=0; ii<block; ii++) {
            if(ii>0)
                printf("    ");
            for(jj=0; jj<block; jj++)
                printf(" %.3e", mat->value[(k*block + ii)*block + jj]);
            printf("\n");
        }
    }
#endif
}

void blockFree(SparseMatrix *mat) {
#ifdef USE_MKL
    mkl_free(mat->value);
    mkl_free(mat->i);
    mkl_free(mat->j);
#else
    free(mat->value);
    free(mat->i);
    free(mat->j);
#endif
}

#elif defined(USE_MKL)
void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName, int wrank);

void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName, int wrank) {
    mat_int k;
    int nread;
    sparse_status_t err;
    sparse_matrix_t coordMatrix;
    mat->i = mkl_malloc(mat->nonzeros * sizeof(*(mat->i)), MALLOC_ALIGN);
    mat->j = mkl_malloc(mat->nonzeros * sizeof(*(mat->j)), MALLOC_ALIGN);
    mat->value = mkl_malloc(mat->nonzeros * sizeof(*(mat->value)), MALLOC_ALIGN);
    for(k=0; k<mat->nonzeros; k++){
        nread = fscanf(fp, "%d%d%le", mat->i+k, mat->j+k, mat->value+k);
        if(nread < 3) {
            fprintf(stderr, "%i:  Error reading %s, element %i\n",
                    wrank, fileName, k);
            break;
        }
    }
    /* Create MKL data structure */
    if((err=mkl_sparse_d_create_coo(&coordMatrix, SPARSE_INDEX_BASE_ONE,
              mat->rows, mat->columns, mat->nonzeros,
                               mat->i, mat->j, mat->value)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "%i:  mkl_sparse_d_create_coo failed %i\n",
                wrank, err);
        exit(55);
    }
    if((err = mkl_sparse_convert_bsr(coordMatrix, mat->blockSize,
                            SPARSE_LAYOUT_COLUMN_MAJOR,
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            &mat->a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "%i:  failed %i\n", wrank, err);
        exit(56);
    }
#if 0  // Unsupported error
    const int expected_calls = 100*1000;
    if((err = mkl_sparse_set_dotmv_hint(mat->a,
            SPARSE_OPERATION_NON_TRANSPOSE, mat->descr, expected_calls)) !=
       SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "%i:  mkl_sparse_set_dotmv_hint failed: %i %i\n",
                wrank, err, SPARSE_STATUS_NOT_SUPPORTED);
        exit(58);
    }
#endif
    if((err = mkl_sparse_optimize(mat->a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "%i:  mkl_sparse_optimize failed %i\n", wrank, err);
        exit(59);
    }
    mkl_free(mat->i); mkl_free(mat->j); mkl_free(mat->value);
    if((err = mkl_sparse_destroy(coordMatrix)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "%i:  failed %i\n", wrank, err);
        exit(60);
    }
 }
#endif

int main(int argc, char **argv){
    char *fdup, *dataPath, *hessFile, *gradFile, *gaugeFile;
    char *options;
    cJSON *jopts, *tmp;
    int i, n, blockSize, gaugeDimension, wsize, wrank;
    mat_int nLargeShifts;
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
    long int twall = 0, tcpu = 0;
    clock_t t1, tt1;
    time_t t2, tt2, tf;

    t1 = clock(); tt1 = t1;
    time(&t2); tt2 = t2;


    /* Initialize MPI */
#ifdef USE_MPI
    if(MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "%i:  Failed to initialize MPI.", wrank);
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
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "blockSize");
    blockSize = cJSON_IsNumber(tmp)?tmp->valueint:1;
    gaugeDimension = cJSON_GetObjectItemCaseSensitive(
                               jopts, "gaugeDimension")->valueint;
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

    /* Read in Hessian Matrix */
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
    assert(hess.rows == abs(n));
    assert(hess.columns == abs(n));
#ifdef USE_BLOCK
    hess.descr = 's';  // Symmetric matrix, with only one triangle stored.
    hess.blockSize = blockSize;
    readtoBlock(fp, hessp, hessFile, wrank);
    sortMatrix(hessp, (chunkSize+1)/2);
#elif defined(USE_MKL)
    hess.descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    hess.descr.mode = SPARSE_FILL_MODE_LOWER;
    hess.descr.diag = SPARSE_DIAG_NON_UNIT;
    hess.blockSize = blockSize;
    readtoBlock(fp, hessp, hessFile, wrank);
#else
    mat_int k;
    hess.descr = 's';  // Symmetric matrix, with only one triangle stored.
    hess.i = malloc(hess.nonzeros * sizeof(*hess.i));
    hess.j = malloc(hess.nonzeros * sizeof(*hess.j));
    hess.value = malloc(hess.nonzeros * sizeof(*hess.value));
    for(k=0; k<hess.nonzeros; k++){
        nread = fscanf(fp, "%u%u%le", hess.i+k, hess.j+k,
                       hess.value+k);
        hess.i[k] -= 1; hess.j[k] -= 1; // switch to zero-based indexing
        if(nread < 3) {
            fprintf(stderr, "%i:  Error reading %s, element %i\n",
                    wrank, hessFile, k);
            break;
        }
    }
    sortMatrix(hessp, (chunkSize+1)/2);
#endif
    fclose(fp);


    mat_int wn;
    rankSanityTest(n/blockSize);

    wn = blockSize*localSize(wrank, wsize, n/blockSize);
    double *grad = malloc(wn * sizeof(double)), *gradp = grad, dummy;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "gradFile");
    gradFile = catStrings(dataPath, tmp->valuestring);
    if(wrank ==0)
        printf("Opening file %s for a vector of length %i\n", gradFile, n);
    fp = fopen(gradFile, "r");
    for(i=0; i<n; i++){
        if(indexRank(i/blockSize, wsize, n/blockSize) == wrank)
            nread = fscanf(fp, "%le", gradp++);
        else
            nread = fscanf(fp, "%le", &dummy);
        if(nread < 1) {
            fprintf(stderr, "%i:  Error reading %s on line %i\n",
                    wrank, gradFile, i);
            break;
        }
    }
    assert(wn == gradp - grad);
    fclose(fp);

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
    assert(gauge.rows == abs(gaugeDimension));
    assert(gauge.columns == abs(n));

#ifdef USE_BLOCK
    gauge.descr = 'g';  // General matrix
    gauge.blockSize = blockSize;
    readtoBlock(fp, gaugep, gaugeFile, wrank);
    sortMatrix(gaugep, chunkSize);
#elif defined(USE_MKL)
    gauge.descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    gauge.descr.mode = SPARSE_FILL_MODE_LOWER;
    gauge.descr.diag = SPARSE_DIAG_NON_UNIT;
    gauge.blockSize = blockSize;
    readtoBlock(fp, gaugep, gaugeFile, wrank);
#else
    gauge.descr = 'g';  // General matrix
    gauge.i = malloc(gauge.nonzeros * sizeof(*gauge.i));
    gauge.j = malloc(gauge.nonzeros * sizeof(*gauge.j));
    gauge.value = malloc(gauge.nonzeros * sizeof(*gauge.value));
    for(k=0; k<gauge.nonzeros; k++){
        nread = fscanf(fp, "%u%u%lf", gauge.i+k, gauge.j+k, gauge.value+k);
        gauge.i[k] -= 1; gauge.j[k] -= 1; // switch to zero-based indexing
        if(nread < 3) {
            fprintf(stderr, "%i:  Error reading %s, element %i\n",
                    wrank, gaugeFile, k);
            break;
        }
    }
    sortMatrix(gaugep, chunkSize);
#endif
    fclose(fp);
    fflush(stdout);
    time(&tf);
    tcpu += clock()-t1;
    twall += tf - t2; 


    /* Solve it!  */

    double *shifts = malloc(n * sizeof(double));
    double *absVals, *vecs;
    int nvals;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "dynamicPartOptions");
    assert(tmp != NULL);
    dynamicInit(gaugep, tmp, mpicomp);
    eigenInit(hessp);
    // Debug print:
#if 0
    testOp(hessp, grad);
#endif
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
    assert(tmp != NULL);
    largeShiftsCheckpoint(grad, tmp, &absVals, &vecs, &nvals, mpicomp);
    cutoffNullspace(n, nvals, jopts, grad, &absVals, &vecs, &nLargeShifts);
    linearInit(hessp, vecs, nLargeShifts);
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "linearSolveOptions");
    assert(tmp != NULL);
    linearSolve(n, grad, tmp, shifts);
    dynamicClose();


    /* output result */

    t1 = clock();
    time(&t2);
    if(wrank ==0)
        printf("Opening output file %s\n", argv[2]);
    fp = fopen(argv[2], "w");
    for(i=0; i<n; i++){
        fprintf(fp, "%.17e\n", shifts[i]);
    }
    fclose(fp);
    time(&tf);
    tcpu += clock()-t1;
    twall += tf - t2;
    printf("%s:  input/output in %.2f sec (%li wall)\n",
           __FILE__, tcpu/(float) CLOCKS_PER_SEC, twall);
    printf("%s:  overall time %.2f sec (%li wall)\n",
           __FILE__, (clock()-tt1)/(float) CLOCKS_PER_SEC, tf-tt2);
    fflush(stdout);


#ifdef USE_BLOCK
    blockFree(hessp);
    blockFree(gaugep);
#elif defined(USE_MKL)
    mkl_sparse_destroy(hess.a); mkl_sparse_destroy(gauge.a);
#else
    free(hess.i); free(hess.j); free(hess.value);
    free(gauge.i); free(gauge.j); free(gauge.value);
#endif
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
