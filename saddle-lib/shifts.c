/*
    Find shifts in link fields for one step 
    of a saddle-point search.

    Example usage:
    ./shifts ../hess-grad-gauge.json ../hess.mtx ../grad.dat ../gauge.mtx ../shifts.dat
    Valgrind debugging example:
    valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./shifts ...

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "shifts.h"
#include "mmio.h"
#ifdef USE_MKL
#include "mkl_spblas.h"
#endif

char *readFile(char *filename) {
    FILE *f = fopen(filename, "rt");
    assert(f);
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buffer = (char *) malloc(length + 1);
    buffer[length] = '\0';
    fread(buffer, 1, length, f);
    fclose(f);
    return buffer;
}

int main(int argc, char **argv){
    char *options;
    cJSON *jopts, *tmp;
    int i, k, n, nc;
    unsigned int nLargeShifts;
    FILE *fp;
#ifdef USE_LIBRSB
    char ib[1000];
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    rsb_flags_t flagsA = RSB_FLAG_NOFLAGS;
    rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
#else
    MM_typecode matcode;
    int j;
#ifdef USE_MKL
    sparse_status_t err;
    const int align = 64;
    int block_size;
    sparse_matrix_t coordMatrix;
#else
    SparseRow *row;
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
#endif

    /* Read JSON file and use options */
    if(argc <5)
        exit(-1);
    printf("Opening file %s", argv[1]);
    options = readFile(argv[1]);
    jopts = cJSON_Parse(options);
    n = cJSON_GetObjectItemCaseSensitive(jopts, "n")->valueint;
    nc = cJSON_GetObjectItemCaseSensitive(jopts, "nc")->valueint;
    printf(" n=%i, nc=%i\n", n, nc);

    /* Read in arrays */
#ifdef USE_LIBRSB
    SparseMatrix *hessp;
    printf("Opening file %s\n", argv[2]);
    hessp = rsb_file_mtx_load(argv[2], flagsA, typecode, &errval);
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
    printf("Opening file %s\n", argv[2]);
    fp = fopen(argv[2], "r"); 
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
    assert(hess.rows == n);
    assert(hess.columns == n);
#ifdef USE_MKL
    hess.row = mkl_malloc(hess.nonzeros * sizeof(int), align);
    hess.column = mkl_malloc(hess.nonzeros * sizeof(int), align);
    hess.value = mkl_malloc(hess.nonzeros * sizeof(double), align);
#else
    hess.data = malloc(hess.nonzeros * sizeof(SparseRow));
#endif
    for(j=0; j<hess.nonzeros; j++){
#ifdef USE_MKL
        k = fscanf(fp, "%d%d%le", hess.row+j, hess.column+j, hess.value+j); 
#else
        row = hess.data+j;
        k = fscanf(fp, "%u%u%le", &(row->i), &(row->j), &(row->value));
        row->i -= 1; row->j -= 1; // switch to zero-based indexing
#endif
        if(k < 3) {
            fprintf(stderr, "Error reading %s, element %i\n", argv[2], j);
            break;
        }
    }
    fclose(fp);
#ifdef USE_MKL
    /* Create MKL data structure */
    if((err=mkl_sparse_d_create_coo(&coordMatrix, SPARSE_INDEX_BASE_ONE,
              hess.rows, hess.columns, hess.nonzeros,
                               hess.row, hess.column, hess.value)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_create_coo failed %i\n", err);
        exit(55);
    }
    block_size = nc*nc - 1;
    if((err = mkl_sparse_convert_bsr(coordMatrix, block_size,
                            SPARSE_LAYOUT_COLUMN_MAJOR, // BvdS: right?
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            &hess.a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "failed %i\n", err);
        exit(56);
    }
    hess.descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    hess.descr.mode = SPARSE_FILL_MODE_LOWER;
    hess.descr.diag = SPARSE_DIAG_NON_UNIT;
#if 0  // Unsupported error
    const int expected_calls = 100*1000;
    if((err = mkl_sparse_set_dotmv_hint(hess.a,
            SPARSE_OPERATION_NON_TRANSPOSE, hess.descr, expected_calls)) !=
       SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_set_dotmv_hint failed: %i %i\n", err,
               SPARSE_STATUS_NOT_SUPPORTED);
        exit(58);
    }
#endif
    if((err = mkl_sparse_optimize(hess.a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_optimize failed %i\n", err);
        exit(59);
    }
    mkl_free(hess.row); mkl_free(hess.column); mkl_free(hess.value);
    if((err = mkl_sparse_destroy(coordMatrix)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "failed %i\n", err);
        exit(60);
    }
#else
    hess.descr = 's';  // Symmetric matrix, with only one triangle stored.
#endif
#endif

    double *grad = malloc(n * sizeof(double));
    printf("Opening file %s for a vector of length %i\n", argv[3], n);
    fp = fopen(argv[3], "r"); 
    for(i=0; i<n; i++){
        k = fscanf(fp, "%le", grad+i);
        if(k < 1) {
            fprintf(stderr, "Error reading %s on line %i\n", argv[3], i);
            break;
        }
    }
    fclose(fp);

#ifdef USE_LIBRSB
    SparseMatrix *gaugep;
    printf("Opening file %s\n", argv[4]);
    gaugep = rsb_file_mtx_load(argv[4], flagsA, typecode, &errval);
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
    printf("Opening file %s\n", argv[4]);
    fp = fopen(argv[4], "r"); 
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner.\n");
        exit(701);
    }
    if(!(mm_is_real(matcode) && mm_is_matrix(matcode) && 
         mm_is_coordinate(matcode) && mm_is_general(matcode))) {
        fprintf(stderr, "Gauge matrix should be real, general, coordinate.\n");
        fprintf(stderr, "Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(702);
    }
    if (mm_read_mtx_crd_size(fp, &(gauge.rows), &(gauge.columns),
                             &(gauge.nonzeros)) !=0) {
        fprintf(stderr, "Cannot read dimensions\n");
        exit(703);
    }
    assert(gauge.columns == n);
#ifdef USE_MKL
    gauge.row = mkl_malloc(gauge.nonzeros * sizeof(int), align);
    gauge.column = mkl_malloc(gauge.nonzeros * sizeof(int), align);
    gauge.value = mkl_malloc(gauge.nonzeros * sizeof(double), align);
#else
    gauge.data = malloc(gauge.nonzeros * sizeof(SparseRow));
#endif
    for(j=0; j<gauge.nonzeros; j++){
#ifdef USE_MKL
        k = fscanf(fp, "%d%d%le", gauge.row+j, gauge.column+j, gauge.value+j); 
#else
        row = gauge.data+j;
        k = fscanf(fp, "%u%u%lf", &(row->i), &(row->j), &(row->value));
        row->i -= 1; row->j -= 1; // switch to zero-based indexing
#endif
        if(k < 3) {
            fprintf(stderr, "Error reading %s, element %i\n", argv[4], j);
            break;
        }
    }
    fclose(fp);
#ifdef USE_MKL
    /* Create MKL data structure */
    if((err = mkl_sparse_d_create_coo(&coordMatrix, SPARSE_INDEX_BASE_ONE,
                gauge.rows, gauge.columns, gauge.nonzeros,
                gauge.row, gauge.column, gauge.value)) !=
            SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_create_coo failed %i\n", err);
        exit(55);
    }
    if((err=mkl_sparse_convert_bsr(coordMatrix, block_size,
             SPARSE_LAYOUT_COLUMN_MAJOR, // BvdS: right?
             SPARSE_OPERATION_NON_TRANSPOSE,
              &gauge.a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_convert_bsr failed %i\n", err);
        exit(56);
    }
    mkl_free(gauge.row); mkl_free(gauge.column); mkl_free(gauge.value);
    if((err=mkl_sparse_destroy(coordMatrix)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_destroy failed %i\n", err);
        exit(57);
    }
    gauge.descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    gauge.descr.mode = SPARSE_FILL_MODE_LOWER;
    gauge.descr.diag = SPARSE_DIAG_NON_UNIT;
#if 0  // Unsupported error
    if((err = mkl_sparse_set_dotmv_hint(gauge.a,
        SPARSE_OPERATION_NON_TRANSPOSE, gauge.descr,
        expected_calls)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_set_dotmv_hint failed %i\n", err);
        exit(58);
    }
#endif
    if((err=mkl_sparse_optimize(gauge.a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_optimize failed %i\n", err);
        exit(58);
    }
#else
    gauge.descr = 'g';  // General matrix
#endif
#endif
    time(&tf);
    tcpu += clock()-t1;
    twall += tf - t2; 

    /* Solve it!  */

    double *shifts = malloc(n * sizeof(double));
    double *absVals, *vecs;
    unsigned int nvals;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "dynamicPartOptions");
    assert(tmp != NULL);
    dynamicInit(gaugep, tmp);
    // Debug print:
#if 0
    testOp(hessp, grad);
#endif
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
    assert(tmp != NULL);
    largeShifts(hessp, grad, tmp, &absVals, &vecs, &nvals);
    cutoffNullspace(n, nvals, jopts, grad, absVals, vecs, &nLargeShifts);
    linearInit(hessp, vecs, nLargeShifts);
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "linearSolveOptions");
    assert(tmp != NULL);
    linearSolve(n, grad, tmp, shifts);
    dynamicClose();


    /* output result */

    t1 = clock();
    time(&t2);
    printf("Opening output file %s\n", argv[5]);
    fp = fopen(argv[5], "w"); 
    for(i=0; i<n; i++){
        fprintf(fp, "%.15e\n", shifts[i]);
    }
    fclose(fp);  
    time(&tf);
    tcpu += clock()-t1;
    twall += tf - t2;
    printf("%s:  input/output in %.2f sec (%li wall)\n",
           __FILE__, tcpu/(float) CLOCKS_PER_SEC, twall);
    printf("%s:  overall time %.2f sec (%li wall)\n",
           __FILE__, (clock()-tt1)/(float) CLOCKS_PER_SEC, tf-tt2);


#ifdef USE_LIBRSB
    printf("rsb_mtx_free: deallocating hess and gauge.\n");
    rsb_mtx_free(hessp);
    rsb_mtx_free(gaugep); 
    rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
#elif defined(USE_MKL)
    mkl_sparse_destroy(hess.a); mkl_sparse_destroy(gauge.a);
#else
    free(hess.data); free(gauge.data);
#endif
    free(absVals); free(vecs); free(shifts); free(grad);
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
