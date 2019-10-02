/*
    Find shifts in link fields for one step 
    of a saddle-point search.

    Example usage:
    ./shifts ../hess-grad-gauge.json ../shifts.dat
    Valgrind debugging example:
    valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ../hess-grad-gauge.json ../shifts.dat

      scp hess-grad-gauge.json hess.mtx gauge.mtx grad.dat bvds@192.168.0.35:lattice/gauge-field-analyze/

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

#ifdef USE_BLOCK
void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName);
int intCmp (const void * a, const void * b);
void blockFree(SparseMatrix *mat);

// See https://stackoverflow.com/questions/10996418
int intCmp (const void * a, const void * b) {
    return (*(size_t *)a < *(size_t *)b) ? -1 :
        (*(size_t *)a > *(size_t *)b);
}
void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName) {
    unsigned int j;
    unsigned int i, k, ii, jj, lastii = 0, blockCols = 0;
    size_t *blockColp;
    const unsigned int block = mat->blockSize,
        maxBlockCol = mat->columns/block;
    double value, *blockRowValue, *blockValuep;
    int *blockColFlag;
    assert(mat->rows%block == 0);
    assert(mat->columns%block == 0);
    mat->blocks = 0;
    blockRowValue = malloc(mat->columns*block*sizeof(double));
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
            k = fscanf(fp, "%u%u%le", &ii, &jj, &value);
            if(k < 3) {
                fprintf(stderr, "Error reading %s, element %i\n", fileName, j);
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
            // Sort the blocks
            qsort(blockColp, blockCols, sizeof(*blockColp), intCmp);
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
void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName);

void readtoBlock(FILE *fp, SparseMatrix *mat, char *fileName) {
    unsigned int j;
    int k;
    sparse_status_t err;
    sparse_matrix_t coordMatrix;
    mat->row = mkl_malloc(mat->nonzeros * sizeof(int), MALLOC_ALIGN);
    mat->column = mkl_malloc(mat->nonzeros * sizeof(int), MALLOC_ALIGN);
    mat->value = mkl_malloc(mat->nonzeros * sizeof(double), MALLOC_ALIGN);
    for(j=0; j<mat->nonzeros; j++){
        k = fscanf(fp, "%d%d%le", mat->row+j, mat->column+j, mat->value+j); 
        if(k < 3) {
            fprintf(stderr, "Error reading %s, element %i\n", fileName, j);
            break;
        }
    }
    /* Create MKL data structure */
    if((err=mkl_sparse_d_create_coo(&coordMatrix, SPARSE_INDEX_BASE_ONE,
              mat->rows, mat->columns, mat->nonzeros,
                               mat->row, mat->column, mat->value)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_create_coo failed %i\n", err);
        exit(55);
    }
    if((err = mkl_sparse_convert_bsr(coordMatrix, mat->blockSize,
                            SPARSE_LAYOUT_COLUMN_MAJOR,
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            &mat->a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "failed %i\n", err);
        exit(56);
    }
#if 0  // Unsupported error
    const int expected_calls = 100*1000;
    if((err = mkl_sparse_set_dotmv_hint(mat->a,
            SPARSE_OPERATION_NON_TRANSPOSE, mat->descr, expected_calls)) !=
       SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_set_dotmv_hint failed: %i %i\n", err,
               SPARSE_STATUS_NOT_SUPPORTED);
        exit(58);
    }
#endif
    if((err = mkl_sparse_optimize(mat->a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_optimize failed %i\n", err);
        exit(59);
    }
    mkl_free(mat->row); mkl_free(mat->column); mkl_free(mat->value);
    if((err = mkl_sparse_destroy(coordMatrix)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "failed %i\n", err);
        exit(60);
    }
 }
#endif

int main(int argc, char **argv){
    char *dataPath, *hessFile, *gradFile, *gaugeFile;
    char *options;
    cJSON *jopts, *tmp;
    int i, k, n, nc, gaugeDimension;
    unsigned int nLargeShifts;
    FILE *fp;
#ifdef USE_LIBRSB
    char ib[1000];
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    rsb_flags_t flagsA = RSB_FLAG_NOFLAGS;
    rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
#else
    MM_typecode matcode;
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
    const int threads = 2;
#endif
#if !defined(USE_BLOCK) && !defined(USE_MKL)
    SparseRow *row;
    unsigned int j;
#endif
#endif
    long int twall = 0, tcpu = 0;
    clock_t t1, tt1;
    time_t t2, tt2, tf;

    t1 = clock(); tt1 = t1;
    time(&t2); tt2 = t2;

#ifdef USE_LIBRSB
    if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_lib_init error 0x%x, exiting\n", errval);
        exit(111);
    }
    /* Doesn't seem to work:  always uses full number of threads,
       However, specifying on the command line works:
             OMP_NUM_THREADS=4 ./shifts ...
       But we see only one CPU fully utilized.  Even a second goes down
       to 50%.  */
    int tn = 2;
    if((errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn))
        != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_lib_set_opt error 0x%x, exiting\n", errval);
        exit(222);
    }
#elif defined(USE_MKL)
    printf("Setting MKL to %i threads\n", threads);
    mkl_set_num_threads_local(threads);
#endif

    /* Read JSON file and use options */
    if(argc <2) {
        fprintf(stderr, "Usage:  ./shifts <parameters.json> <output.dat>\n");
        exit(-1);
    }
    dataPath = dirname(strdup(argv[1]));
    printf("Opening file %s", argv[1]);
    options = readFile(argv[1]);
    jopts = cJSON_Parse(options);
    n = cJSON_GetObjectItemCaseSensitive(jopts, "n")->valueint;
    nc = cJSON_GetObjectItemCaseSensitive(jopts, "nc")->valueint;
    gaugeDimension = cJSON_GetObjectItemCaseSensitive(
                               jopts, "gaugeDimension")->valueint;
    printf(" n=%i, nc=%i, gauge=%i\n", n, nc, gaugeDimension);

    /* Read in Hessian Matrix */
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "hessFile");
    hessFile = catStrings(dataPath, tmp->valuestring);
#ifdef USE_LIBRSB
    SparseMatrix *hessp;
    printf("Opening file %s\n", hessFile);
    hessp = rsb_file_mtx_load(hessFile, flagsA, typecode, &errval);
    if(errval != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_file_mtx_load error 0x%x, exiting\n", errval);
        exit(113);
    }
    /* print out the matrix summary information  */
    rsb_mtx_get_info_str(hessp, "RSB_MIF_MATRIX_INFO__TO__CHAR_P",
                         ib, sizeof(ib));
    printf("%s\n",ib);
#else
    SparseMatrix hess, *hessp = &hess;
    printf("Opening file %s", hessFile);
    if((fp = fopen(hessFile, "r")) == NULL) {
        fprintf(stderr, "Could not open file %s\n", hessFile);
        exit(44);
    }
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner.\n");
        exit(701);
    }
    if(!(mm_is_real(matcode) && mm_is_matrix(matcode) && 
         mm_is_coordinate(matcode))) {
        fprintf(stderr, "Hessian should be real, coordinate.\n");
        fprintf(stderr, "Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(702);
    }
    if (mm_read_mtx_crd_size(fp, &(hess.rows), &(hess.columns),
                             &(hess.nonzeros)) !=0) {
        fprintf(stderr, "Cannot read dimensions\n");
        exit(703);
    }
    printf(" with %i nonzero elements.\n", hess.nonzeros);
    assert(hess.rows == abs(n));
    assert(hess.columns == abs(n));
#ifdef USE_BLOCK
    hess.descr = 's';  // Symmetric matrix, with only one triangle stored.
    hess.blockSize = nc*nc-1;
    readtoBlock(fp, hessp, hessFile);
#elif defined(USE_MKL)
    hess.descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    hess.descr.mode = SPARSE_FILL_MODE_LOWER;
    hess.descr.diag = SPARSE_DIAG_NON_UNIT;
    hess.blockSize = nc*nc-1;
    readtoBlock(fp, hessp, hessFile);
#else
    hess.descr = 's';  // Symmetric matrix, with only one triangle stored.
    hess.data = malloc(hess.nonzeros * sizeof(SparseRow));
    for(j=0; j<hess.nonzeros; j++){
        row = hess.data+j;
        k = fscanf(fp, "%u%u%le", &(row->i), &(row->j), &(row->value));
        row->i -= 1; row->j -= 1; // switch to zero-based indexing
        if(k < 3) {
            fprintf(stderr, "Error reading %s, element %i\n", hessFile, j);
            break;
        }
    }
#endif
    fclose(fp);
#endif

    double *grad = malloc(n * sizeof(double));
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "gradFile");
    gradFile = catStrings(dataPath, tmp->valuestring);
    printf("Opening file %s for a vector of length %i\n", gradFile, n);
    fp = fopen(gradFile, "r");
    for(i=0; i<n; i++){
        k = fscanf(fp, "%le", grad+i);
        if(k < 1) {
            fprintf(stderr, "Error reading %s on line %i\n", gradFile, i);
            break;
        }
    }
    fclose(fp);

    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "gaugeFile");
    gaugeFile = catStrings(dataPath, tmp->valuestring);
#ifdef USE_LIBRSB
    SparseMatrix *gaugep;
    printf("Opening file %s\n", gaugeFile);
    gaugep = rsb_file_mtx_load(gaugeFile, flagsA, typecode, &errval);
    if(errval != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_file_mtx_load error 0x%x, exiting\n", errval);
        exit(113);
    }
    /* print out the matrix summary information  */
    rsb_mtx_get_info_str(gaugep, "RSB_MIF_MATRIX_INFO__TO__CHAR_P",
                         ib, sizeof(ib));
    printf("%s\n",ib);
#else
    SparseMatrix gauge, *gaugep = &gauge;
    printf("Opening file %s", gaugeFile);
    fp = fopen(gaugeFile, "r");
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner.\n");
        exit(701);
    }
    if(!(mm_is_real(matcode) && mm_is_matrix(matcode) && 
         mm_is_coordinate(matcode) && mm_is_general(matcode))) {
        fprintf(stderr, "Gauge matrix should be real, "
                         "general, coordinate.\n");
        fprintf(stderr, "Market Market type: [%s]\n",
                mm_typecode_to_str(matcode));
        exit(702);
    }
    if (mm_read_mtx_crd_size(fp, &(gauge.rows), &(gauge.columns),
                             &(gauge.nonzeros)) !=0) {
        fprintf(stderr, "Cannot read dimensions\n");
        exit(703);
    }
    printf(" with %i nonzero elements.\n", gauge.nonzeros);
    assert(gauge.rows == abs(gaugeDimension));
    assert(gauge.columns == abs(n));

#ifdef USE_BLOCK
    gauge.descr = 'g';  // General matrix
    gauge.blockSize = nc*nc-1;
    readtoBlock(fp, gaugep, gaugeFile);
#elif defined(USE_MKL)
    gauge.descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    gauge.descr.mode = SPARSE_FILL_MODE_LOWER;
    gauge.descr.diag = SPARSE_DIAG_NON_UNIT;
    gauge.blockSize = nc*nc-1;
    readtoBlock(fp, gaugep, gaugeFile);
#else
    gauge.descr = 'g';  // General matrix
    gauge.data = malloc(gauge.nonzeros * sizeof(SparseRow));
    for(j=0; j<gauge.nonzeros; j++){
        row = gauge.data+j;
        k = fscanf(fp, "%u%u%lf", &(row->i), &(row->j), &(row->value));
        row->i -= 1; row->j -= 1; // switch to zero-based indexing
        if(k < 3) {
            fprintf(stderr, "Error reading %s, element %i\n", gaugeFile, j);
            break;
        }
    }
#endif
    fclose(fp);
#endif
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
    dynamicInit(gaugep, tmp);
    eigenInit(hessp);
    // Debug print:
#if 0
    testOp(hessp, grad);
#endif
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
    assert(tmp != NULL);
    largeShiftsCheckpoint(grad, tmp, &absVals, &vecs, &nvals);
    cutoffNullspace(n, nvals, jopts, grad, &absVals, &vecs, &nLargeShifts);
    linearInit(hessp, vecs, nLargeShifts);
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "linearSolveOptions");
    assert(tmp != NULL);
    linearSolve(n, grad, tmp, shifts);
    dynamicClose();


    /* output result */

    t1 = clock();
    time(&t2);
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


#ifdef USE_LIBRSB
    printf("rsb_mtx_free: deallocating hess and gauge.\n");
    rsb_mtx_free(hessp);
    rsb_mtx_free(gaugep); 
    rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
#elif defined(USE_BLOCK)
    blockFree(hessp);
    blockFree(gaugep);
#elif defined(USE_MKL)
    mkl_sparse_destroy(hess.a); mkl_sparse_destroy(gauge.a);
#else
    free(hess.data); free(gauge.data);
#endif
    free(dataPath);
    free(absVals); free(vecs); free(shifts); free(grad);
    free(hessFile); free(gradFile); free(gaugeFile);
    cJSON_Delete(jopts);
    free(options);
    return 0;
}

#ifdef USE_LIBRSB
int rows(SparseMatrix *matrix){
    int value;
    rsb_mtx_get_info(matrix,
                     RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &value);
    return value;
}
int columns(SparseMatrix *matrix){
    int value;
    rsb_mtx_get_info(matrix,
                     RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &value);
    return value;
}
int nonzeros(SparseMatrix *matrix){
    int value;
    rsb_mtx_get_info(matrix,
                     RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &value);
    return value;
}
#endif
