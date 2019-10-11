/*
  Sort block-sparse coordinate format matrix into
  further blocks, "chunks," to minimize cache misses.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "shifts.h"

void makeChunkNumber(const SparseMatrix *mat,
                     const mat_int chunk);
int indexCmp (const void * a, const void * b);

mat_int *chunkNumber;

/* Calculate index for each block of the matrix,
   so that "chunks" are continguous.
   
   Note that the number of block rows mat->rows/mat->blockSize 
   need not be a multiple of chunk.  In that case,
   there will be holes in the numbering, but the numbering will
   still be correctly ordered.
*/
void makeChunkNumber(const SparseMatrix *mat,
                         const mat_int chunk) {
    mat_int i, j, k;
#ifdef USE_BLOCK
    const mat_int b1 = mat->blockSize, n=mat->blocks;
#else
    const mat_int b1 = 1, n = mat->nonzeros;
#endif
    assert(chunk>0);
    const mat_int columnChunks =
        (mat->columns/b1 + chunk - 1)/chunk;

    chunkNumber = malloc(n * sizeof(*chunkNumber));

    for(k=0; k<n; k++) {
        i = (mat->i)[k]/b1;
        j = (mat->j)[k]/b1;
        chunkNumber[k] = (((i/chunk)*columnChunks +
            j/chunk)*chunk + i%chunk)*chunk + j%chunk;
    }
}

/* See https://stackoverflow.com/questions/10996418
   https://stackoverflow.com/questions/31229657 */
int indexCmp (const void * a, const void * b) {
    const mat_int ia = *(mat_int *)a,
                       ib = *(mat_int *)b;
    return (chunkNumber[ia] > chunkNumber[ib]) -
        (chunkNumber[ib] > chunkNumber[ia]);
}


/*
  Sort matrix blocks so into chunk-size blocks.
  For the right value of chunk, cache misses
  should be minimized
*/
void sortMatrix(SparseMatrix *mat, const mat_int chunk) {
    mat_int k, k1, k2, *index;
    mat_int holdi, holdj;
#ifdef USE_BLOCK
    const mat_int n = mat->blocks,
        nn = mat->blockSize * mat->blockSize;
#else
    const mat_int n = mat->nonzeros,
        nn = 1;
#endif
    double *holdValue = malloc(nn*sizeof(*holdValue));
    index = malloc(n*sizeof(*index));
    for(k=0; k<n; k++)
        index[k] = k;

    // Create index by sorting with respect to chunkNumber
    makeChunkNumber(mat, chunk);
    qsort(index, n, sizeof(*index), indexCmp);
#if 0
    // Print mathematica format coordinates
printf("preCoords = {");
    for(k=0; k<n; k++) {
        if(k>0)
            printf(", ");
        if(k%10==-1)
            printf("\n");
printf("{%i, %i, %i}", mat->i[k], mat->j[k],
           chunkNumber[k]);
    }
    printf("};\n");
#endif
    free(chunkNumber);

    /* Reorder matrix blocks in-place based on index
       See https://stackoverflow.com/questions/7365814 */
    //  i_dst_first -> k
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
        memcpy(holdValue, mat->value+k2*nn,
                   nn*sizeof(*(mat->value)));

        // Follow permutation cycle
        do {
            mat->i[k2] = mat->i[k1];
            mat->j[k2] = mat->j[k1];
            memcpy(mat->value+k2*nn, mat->value+k1*nn,
                       nn*sizeof(*(mat->value)));
            index[k2] = k2;

            k2 = k1;
            k1 = index[k1];
            assert(k1 != k2);
        } while(k1 != k);

        mat->i[k2] = holdi;
        mat->j[k2] = holdj;
        memcpy(mat->value+k2*nn, holdValue,
                   nn*sizeof(*(mat->value)));
        index[k2] = k2;
    }

#if 0
    // Print mathematica format coordinates
printf("postCoords = {");
    for(k=0; k<n; k++) {
        if(k>0)
            printf(", ");
        if(k%10==-1)
            printf("\n");
printf("{%i, %i, %i}", mat->i[k], mat->j[k], k);
    }
    printf("};\n");
#endif

    free(holdValue);
    free(index);
}
