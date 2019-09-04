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
    long int twall = 0, tcpu = 0;
    clock_t t1, tt1;
    time_t t2, tt2, tf;

    t1 = clock(); tt1 = t1;
    time(&t2); tt2 = t2;

    /* Read JSON file and use options */
    if(argc <5)
        exit(-1);
    printf("Opening file %s\n", argv[1]);
    options = readFile(argv[1]);
    jopts = cJSON_Parse(options);
    n = cJSON_GetObjectItemCaseSensitive(jopts, "n")->valueint;

    /* Read in arrays */
    unsigned int hessElements = cJSON_GetObjectItemCaseSensitive(
             jopts, "hessElements")->valueint;
    SparseRow *hess = malloc(hessElements * sizeof(SparseRow));
    printf("Opening file %s for %i elements\n", argv[2], hessElements);
    fp = fopen(argv[2], "r"); 
    for(i=0; i<hessElements; i++){
        k = fscanf(fp, "%u%u%le", &(hess[i].i), &(hess[i].j), &(hess[i].value));
        if(k < 3) {
            printf("Error reading %s on line %i\n", argv[2], i);
            break;
        }
    }
    fclose(fp);
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
    unsigned int gaugeElements = cJSON_GetObjectItemCaseSensitive(
             jopts, "gaugeElements")->valueint;
    unsigned int gaugeDimension = cJSON_GetObjectItemCaseSensitive(
             jopts, "gaugeDimension")->valueint;
    SparseRow *gauge = malloc(gaugeElements * sizeof(SparseRow));
    printf("Opening file %s for %i elements\n", argv[4], gaugeElements);
    fp = fopen(argv[4], "r"); 
    for(i=0; i<gaugeElements; i++){
        k = fscanf(fp, "%u%u%lf", &(gauge[i].i), &(gauge[i].j), &(gauge[i].value));
        if(k < 3) {
            printf("Error reading %s on line %i\n", argv[4], i);
            break;
        }
    }
    fclose(fp);
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
    dynamicInit(n, gauge, gaugeDimension, gaugeElements, tmp);
    hessInit(n, hess, hessElements);
    // Debug print:
#if 0
    testOp(n, grad);
#endif
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
    assert(tmp != NULL);
    /* This won' twork until we introduce reorthogonalization against
       gauge shifts in TrLAN */
    largeShifts(n, grad, tmp, &vals, &vecs, &nvals);
    cutoffNullspace(n, nvals, jopts, grad, vals, vecs, &nLargeShifts);
    linearInit(n, hess, hessElements, vecs, nLargeShifts);
    tmp = cJSON_GetObjectItemCaseSensitive(jopts, "linearSolveOptions");
    assert(tmp != NULL);
    linearSolve(n, grad, tmp, shifts);
    printDynamicStats();
  
    /* output result */
    t1 = clock();
    time(&t2);
    printf("Opening file %s\n", argv[5]);
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

    free(hess); free(grad); free(gauge);
    free(vals); free(vecs); free(shifts);
    cJSON_Delete(jopts);
    free(options);
    return 0;
}
