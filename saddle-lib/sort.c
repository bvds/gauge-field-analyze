/*
  Sort block-sparse coordinate format matrix into
  further blocks, "chunks," to minimize cache misses.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "shifts.h"

void doSortMatrix(SparseMatrix *mat, int debug);
int indexCmp (const void * a, const void * b);
void printCoordinates(SparseMatrix *mat, char* variableName, int rankFlag);

// Type must accommodate n^2
typedef unsigned long mat_index;
mat_index *elementRank;

/* Calculate index for each block of the matrix,
   so that "chunks" are continguous.
   
   Note that the number of block rows mat->rows/mat->blockSize 
   need not be a multiple of chunk.  In that case,
   there will be holes in the numbering, but the numbering will
   still be correctly ordered.
*/
void sortMatrixChunks(SparseMatrix *mat,
                         const mat_int chunkIn) {
    mat_index i, j;
    mat_int k;
    const mat_index chunk = chunkIn;
    const mat_int b1 = mat->blockSize, n = mat->blocks;
    assert(chunk>0);
    const mat_index columnChunks =
        (mat->columns/b1 + chunk - 1)/chunk;
    elementRank = malloc(n * sizeof(*elementRank));

    for(k=0; k<n; k++) {
        i = (mat->i)[k]/b1;
        j = (mat->j)[k]/b1;
        elementRank[k] = (((i/chunk)*columnChunks +
            j/chunk)*chunk + i%chunk)*chunk + j%chunk;
    }

    doSortMatrix(mat, 0);
    free(elementRank);
}


/*
  Sort matrix into blocks based on MPI processes.
  Shift columns so that blocks associated with the
  same process come first.

  This assumes that the rows have been filtered
  to the local process and the row indices have been
  shifted.

  For a single process, this should just sort blocks
  by row and column.
*/
void sortMatrixLocal(SparseMatrix *mat, int wrank, int wsize, int debug) {
    mat_index i, j, j0;
    mat_int k;
    int jrank;
    const mat_index cols = mat->columns;
    const mat_int n = mat->blocks;
    const mat_index rowPart = mat->rows/mat->rowParts,
        colPart = mat->columns/mat->colParts;
    const mat_index localRows = rowPart*localSize(wrank, wsize, mat->rowParts);
    const mat_index localj0 = colPart*rankIndex(wrank, wsize, mat->colParts);

    elementRank = malloc(n * sizeof(*elementRank));

    for(k=0; k<n; k++) {
        i = (mat->i)[k];
        j = (mat->j)[k];
        assert(i < localRows);
        jrank = indexRank(j/colPart, wsize, mat->colParts);
        j0 = colPart*rankIndex(jrank, wsize, mat->colParts);
        /* C standard for mod of negative numbers is screwy.
           add cols to avoid issue. */
        elementRank[k] = localRows*((cols + j0 - localj0)%cols) +
            i*colPart*localSize(jrank, wsize, mat->colParts) +
            j - j0;
        assert(elementRank[k] < localRows*cols);
    }

    doSortMatrix(mat, debug);
    free(elementRank);
}

/* See https://stackoverflow.com/questions/10996418
   https://stackoverflow.com/questions/31229657 */
int indexCmp (const void * a, const void * b) {
    const mat_int ia = *(mat_int *)a,
        ib = *(mat_int *)b;
    return (elementRank[ia] > elementRank[ib]) -
        (elementRank[ib] > elementRank[ia]);
}


/*
  Sort matrix blocks so based on a given elementRank.
*/
void doSortMatrix(SparseMatrix *mat, int debug) {
    mat_int k, k1, k2, *index;
    mat_int holdi, holdj;
    const mat_int n = mat->blocks,
        b2 = mat->blockSize * mat->blockSize;
    double *holdValue = malloc(b2*sizeof(*holdValue));

    if(debug>1)
        printCoordinates(mat, "preCoords", 1); // elementRank label
    // Create index by sorting with respect to elementRank
    index = malloc(n*sizeof(*index));
    for(k=0; k<n; k++)
        index[k] = k;
    qsort(index, n, sizeof(*index), indexCmp);

    /* Reorder matrix blocks in-place based on index
       See https://stackoverflow.com/questions/7365814 */
    // i_dst_first -> k
    // i_src -> k1
    // i_dst -> k2
    for(k=0; k<n; k++) {
        k1 = index[k];
        assert(k1 < n);
        if(k1 == k)
            continue; // already sorted
        k2 = k;
        holdi = mat->i[k2];
        holdj = mat->j[k2];
        memcpy(holdValue, mat->value+k2*b2,
                   b2*sizeof(*(mat->value)));

        // Follow permutation cycle
        do {
            mat->i[k2] = mat->i[k1];
            mat->j[k2] = mat->j[k1];
            memcpy(mat->value+k2*b2, mat->value+k1*b2,
                       b2*sizeof(*(mat->value)));
            index[k2] = k2;

            k2 = k1;
            k1 = index[k1];
            assert(k1 != k2);
        } while(k1 != k);

        mat->i[k2] = holdi;
        mat->j[k2] = holdj;
        memcpy(mat->value+k2*b2, holdValue,
                   b2*sizeof(*(mat->value)));
        index[k2] = k2;
    }

    if(debug>0)
        printCoordinates(mat, "postCoords", 0); // actual order

    free(holdValue);
    free(index);
}


/*
   debug print
   print coordinates in Mathematica format.
*/
void printCoordinates(SparseMatrix *mat, char* variableName, int rankFlag) {
    mat_int k;
    // Print mathematica format coordinates
    printf("blockSize=%i;\n", mat->blockSize);
    printf("%s = {", variableName);
    for(k=0; k<mat->blocks; k++) {
        if(k>0)
            printf(", ");
        if(k%10==9)
            printf("\n");
        printf("{%i, %i, %li}", mat->i[k], mat->j[k],
               rankFlag?elementRank[k]:k);
    }
    printf("};\n");
}
