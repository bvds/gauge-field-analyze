#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "shifts.h"
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


#ifdef USE_BLOCK
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

/* in and out must be distinct */
void matrixVector(const SparseMatrix *a,
                  const doublereal *in, doublereal *out) {
#ifdef USE_BLOCK
    size_t k;
    double *matp = a->value;
    const integer n = a->blockSize, n2=n*n;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T', normal='N';

    memset(out, 0, a->rows * sizeof(double));
    for(k=0; k<a->blocks; k++, matp+=n2) {
        DGEMV(&trans, &n, &n, &one,
              matp, &n, in + a->j[k], &inc, &one,
              out + a->i[k], &inc);

        if(a->descr == 's' && a->i[k] != a->j[k]) {
            DGEMV(&normal, &n, &n, &one,
                  matp, &n, in + a->i[k], &inc, &one,
                  out + a->j[k], &inc);
        }
    }
#elif defined(USE_MKL)
    sparse_status_t err;
    const double alpha = 1.0, beta=0.0;
    double d;

    if((err=mkl_sparse_d_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
              alpha, a->a, a->descr, in, beta, out, &d)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_dotmv failed %i\n", err);
        exit(69);
    }
#else
    mat_int i, j, k;
    double value;
    memset(out, 0, a->rows * sizeof(double));
    for(k=0; k<a->nonzeros; k++) {
        i = a->i[k];
        j = a->j[k];
        value = a->value[k];
        out[i] += value * in[j];
        if(a->descr == 's' && j != i) {
            out[j] += value * in[i];
        }
    }
#endif
}

void vectorMatrix(const SparseMatrix *a,
                  const doublereal *in, doublereal *out) {
#ifdef USE_BLOCK
    size_t k;
    double *matp = a->value;
    const integer n = a->blockSize, n2 = n*n;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T', normal='N';

    memset(out, 0, a->columns * sizeof(*out));
    for(k=0; k<a->blocks; k++, matp+=n2) {
        DGEMV(&normal, &n, &n, &one,
              matp, &n, in + a->i[k], &inc, &one,
              out + a->j[k], &inc);

        if(a->descr == 's' && a->i[k] != a->j[k]) {
            DGEMV(&trans, &n, &n, &one,
                  matp, &n, in + a->j[k], &inc, &one,
                  out + a->i[k], &inc);
        }
    }
#elif defined(USE_MKL)
    sparse_status_t err;
    const double alpha = 1.0, beta=0.0;
    double d;

    if((err=mkl_sparse_d_dotmv(SPARSE_OPERATION_TRANSPOSE,
              alpha, a->a, a->descr, in, beta, out, &d)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_dotmv failed %i\n", err);
        exit(68);
    }
#else
    mat_int i, j, k;
    double value;
    memset(out, 0, a->columns * sizeof(*out));
    for(k=0; k<a->nonzeros; k++) {
        i = a->i[k];
        j = a->j[k];
        value = a->value[k];
        out[j] += value * in[i];
        if(a->descr == 's' && i != j) {
            out[i] += value * in[j];
        }
    }
#endif
}
