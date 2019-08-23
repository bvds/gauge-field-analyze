/*
    Find shifts in link fields for one step 
    in saddle-point search.

    Example usage:
    ./shifts ../hess-grad-gauge.json ../hess.dat ../grad.dat ../gauge.dat ../shifts.dat

    Use https://github.com/DaveGamble/cJSON
    installed locally.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <cjson/cJSON.h>

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

struct SparseRow {
  unsigned int i;
  unsigned int j;
  double value;
};

int main(int argc, char **argv){
  char *options;
  cJSON *jopts, *tmp;
  unsigned int n, nc;
  double rtol1;
  int i;
  FILE *fp;

  /* Read JSON file and use options */
  printf("Opening file %s\n", argv[1]);
  options = readFile(argv[1]);
  jopts = cJSON_Parse(options);
  n = cJSON_GetObjectItemCaseSensitive(jopts, "n")->valueint;
  nc = cJSON_GetObjectItemCaseSensitive(jopts, "nc")->valueint;
  printf("n=%i Nc = %i\n", n, nc);
  /* Example nested get */
  tmp = cJSON_GetObjectItemCaseSensitive(jopts, "largeShiftOptions");
  rtol1 = cJSON_GetObjectItemCaseSensitive(tmp, "rTolerance")->valuedouble;
  printf("rtol1=%e\n", rtol1);

  /* Read in arrays */
  int hessElements = cJSON_GetObjectItemCaseSensitive(jopts,
				"hessElements")->valueint;
  struct SparseRow *hess = malloc(hessElements * sizeof(struct SparseRow));
  printf("Opening file %s\n", argv[2]);
  fp = fopen(argv[2], "r"); 
  for(i=0; i<hessElements; i++){
    fscanf(fp, "%i %i %le", &(hess[i].i), &(hess[i].j), &(hess[i].value));
  }
  fclose(fp);
  double *grad = malloc(n * sizeof(double));
  printf("Opening file %s\n", argv[3]);
  fp = fopen(argv[3], "r"); 
  for(i=0; i<n; i++){
    fscanf(fp, "%le", grad+i);
  }
  fclose(fp);
  int gaugeElements = cJSON_GetObjectItemCaseSensitive(jopts,
				"gaugeElements")->valueint;
  struct SparseRow *gauge = malloc(gaugeElements * sizeof(struct SparseRow));
  printf("Opening file %s\n", argv[4]);
  fp = fopen(argv[4], "r"); 
  for(i=0; i<gaugeElements; i++){
    fscanf(fp, "%i %i %le", &(gauge[i].i), &(gauge[i].j), &(gauge[i].value));
  }
  fclose(fp);

  /* Solve it!  */

  double *shifts = malloc(n * sizeof(double));
  for(i=0; i<n; i++){
    shifts[i] = 0.5;
  }
  
  /* output result */
  printf("Opening file %s\n", argv[5]);
  fp = fopen(argv[5], "w"); 
  for(i=0; i<n; i++){
    fprintf(fp, "%.15e\n", shifts[i]);
  }
  fclose(fp);  

  free(hess); free(grad); free(gauge);
  cJSON_Delete(jopts);
  free(options);
  return 0;
}
