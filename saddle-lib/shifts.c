/*
    Find shifts in link fields for one step 
    in saddle-point search.

    Example usage:
    ./shifts ../hess-grad-gauge.json ../hess.dat ../grad.dat ../gauge.dat ../shifts.dat
    Valgrind debugging example:
    valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./shifts ...

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "shifts.h"

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
    unsigned int n;
    unsigned int i;
    unsigned int nLargeShifts;
    int k;
    FILE *fp;
#ifdef USE_LIBRSB
    char ib[1000];
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    rsb_flags_t flagsA = RSB_FLAG_NOFLAGS;
    rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
#else
    SparseRow *row;
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
#endif

    /* Read JSON file and use options */
    if(argc <5)
        exit(-1);
    printf("Opening file %s\n", argv[1]);
    options = readFile(argv[1]);
    jopts = cJSON_Parse(options);
    n = cJSON_GetObjectItemCaseSensitive(jopts, "n")->valueint;

    /* Read in arrays */
    SparseMatrix hess;
#ifdef USE_LIBRSB
    printf("Opening file %s\n", argv[2]);
    hess = rsb_file_mtx_load(argv[2], flagsA, typecode, &errval);
    if(errval != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_file_mtx_load error 0x%x, exiting\n", errval);
        exit(113);
    }
    /* print out the matrix summary information  */
    rsb_mtx_get_info_str(hess, "RSB_MIF_MATRIX_INFO__TO__CHAR_P",
                         ib, sizeof(ib));
    printf("%s\n",ib);
#else
    hess.nonzeros = cJSON_GetObjectItemCaseSensitive(
             jopts, "hessElements")->valueint;
    hess.rows = n;
    hess.columns = n;
    hess.data = malloc(hess.nonzeros * sizeof(SparseRow));
    printf("Opening file %s for %i elements\n", argv[2], hess.nonzeros);
    fp = fopen(argv[2], "r"); 
    for(i=0; i<hess.nonzeros; i++){
        row = hess.data+i;
        k = fscanf(fp, "%u%u%le", &(row->i), &(row->j), &(row->value));
        if(k < 3) {
            printf("Error reading %s on line %i\n", argv[2], i);
            break;
        }
    }
    fclose(fp);
#endif

    double *grad = malloc(n * sizeof(double));
    printf("Opening file %s for %i elements\n", argv[3], n);
    fp = fopen(argv[3], "r"); 
    for(i=0; i<n; i++){
        k = fscanf(fp, "%le", grad+i);
        if(k < 1) {
            printf("Error reading %s on line %i\n", argv[3], i);
            break;
        }
    }
    fclose(fp);
    SparseMatrix gauge;
#ifdef USE_LIBRSB
    printf("Opening file %s\n", argv[4]);
    gauge = rsb_file_mtx_load(argv[4], flagsA, typecode, &errval);
    if(errval != RSB_ERR_NO_ERROR) {
        fprintf(stderr, "rsb_file_mtx_load error 0x%x, exiting\n", errval);
        exit(113);
    }
    /* print out the matrix summary information  */
    rsb_mtx_get_info_str(gauge, "RSB_MIF_MATRIX_INFO__TO__CHAR_P",
                         ib, sizeof(ib));
    printf("%s\n",ib);
#else
    gauge.nonzeros = cJSON_GetObjectItemCaseSensitive(
             jopts, "gaugeElements")->valueint;
    gauge.rows = cJSON_GetObjectItemCaseSensitive(
             jopts, "gaugeDimension")->valueint;
    gauge.columns = n;
    gauge.data = malloc(gauge.nonzeros * sizeof(SparseRow));
    printf("Opening file %s for %i elements\n", argv[4], gauge.nonzeros);
    fp = fopen(argv[4], "r"); 
    for(i=0; i<gauge.nonzeros; i++){
        row = gauge.data+i;
        k = fscanf(fp, "%u%u%lf", &(row->i), &(row->j), &(row->value));
        if(k < 3) {
            printf("Error reading %s on line %i\n", argv[4], i);
            break;
        }
    }
    fclose(fp);
#endif
    time(&tf);
    tcpu += clock()-t1;
    twall += tf - t2; 

    /* Solve it!  */

    /* Find lowest eigenpairs, but do nothing with them. */
    double *shifts = malloc(n * sizeof(double));
    double *vals, *vecs;
    unsigned int nvals;
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "dynamicPartOptions");
    assert(tmp != NULL);
    dynamicInit(&gauge, tmp);
    // Debug print:
#if 0
    testOp(&hess, grad);
#endif
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
    assert(tmp != NULL);
    /* This won' twork until we introduce reorthogonalization against
       gauge shifts in TrLAN */
    largeShifts(&hess, grad, tmp, &vals, &vecs, &nvals);
    cutoffNullspace(n, nvals, jopts, grad, vals, vecs, &nLargeShifts);
    linearInit(&hess, vecs, nLargeShifts);
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
    rsb_mtx_free(hess);
    rsb_mtx_free(gauge); 
    rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
#else
    free(hess.data); free(gauge.data);
#endif
    free(vals); free(vecs); free(shifts); free(grad);
    cJSON_Delete(jopts);
    free(options);
    return 0;
}
