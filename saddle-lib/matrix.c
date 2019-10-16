#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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


void sparseMatrixRead(FILE *fp, SparseMatrix *mat, char *fileName,
                      int blockSize, _MPI_Comm mpicom) {
    mat_int k;
    int nread;
#ifdef USE_MPI
    int wrank, wsize;

    mat->mpicom = mpicom;
    MPI_Comm_rank(mpicom, &wrank);
    MPI_Comm_size(mpicom, &wsize);
#else
    int wrank = 0;
    assert(mpicom == NULL);
#endif
#ifdef USE_BLOCK
    mat_int i, j, ii, jj, lastii = 0, blockCols = 0;
    mat_int *blockColp;
    const mat_int maxBlockCol = mat->columns/blockSize;
    double value, *blockRowValue, *blockValuep;
    int *blockColFlag;
    assert(mat->rows%blockSize == 0);
    assert(mat->columns%blockSize == 0);

    mat->blockSize = blockSize;
    mat->blocks = 0;
    blockRowValue = malloc(mat->columns*blockSize*sizeof(*blockRowValue));
    blockColFlag = malloc(maxBlockCol*sizeof(*blockColFlag));
    blockColp = malloc(maxBlockCol*sizeof(*blockColp));
    mat->value = MALLOC(0);
    mat->i = MALLOC(0);
    mat->j = MALLOC(0);
    memset(blockRowValue, 0, mat->columns*blockSize*sizeof(double));
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
        if(j==mat->nonzeros || (lastii<ii && ii%blockSize == 0)) {
            // Retire block row and reset for a new one
            mat->value = REALLOC(mat->value,
                     (mat->blocks + blockCols)*blockSize*blockSize*sizeof(double));
            mat->i = REALLOC(mat->i,
                     (mat->blocks + blockCols)*sizeof(*mat->i));
            mat->j = REALLOC(mat->j,
                     (mat->blocks + blockCols)*sizeof(*mat->j));
            for(i=0; i<blockCols; i++) {
                k = blockColp[i];
#if 0
                printf("  %i %i mat->blocks=%i blockCols=%i\n",
                       j, k, mat->blocks, blockCols);
#endif
                blockValuep = blockRowValue+k*blockSize*blockSize;
                memcpy(mat->value+blockSize*blockSize*mat->blocks,
                       blockValuep, blockSize*blockSize*sizeof(double));
                memset(blockValuep, 0, blockSize*blockSize*sizeof(double));
                mat->i[mat->blocks] = (lastii/blockSize)*blockSize;
                mat->j[mat->blocks] = k*blockSize;
                mat->blocks += 1;
                blockColFlag[k] = 0;
            }
            blockCols = 0;
        }
        if(j==mat->nonzeros)
            break;
        lastii = ii;
        // Mark this block
        k = jj/blockSize;
        if(!blockColFlag[k]) {
            blockColp[blockCols] = k;
            blockCols += 1;
            blockColFlag[k] = 1;
        }
        // Add value to this block
        blockRowValue[(k*blockSize + ii%blockSize)*blockSize + jj%blockSize] = value;
        /* fill diagonal blocks of symmetric matrix */
        if(mat->descr == 's' && ii/blockSize == k && ii!=jj) {
            blockRowValue[(k*blockSize + jj%blockSize)*blockSize + ii%blockSize] = value;
        }
    }
    free(blockRowValue);
    free(blockColFlag);
    free(blockColp);
#if 0
    printf("Print matrix for file %s\n", fileName);
    for(k=0; k<mat->blocks; k++) {
        printf("%i %i  ", mat->i[k] + 1, mat->j[k] + 1);
        for(ii=0; ii<blockSize; ii++) {
            if(ii>0)
                printf("    ");
            for(jj=0; jj<blockSize; jj++)
                printf(" %.3e", mat->value[(k*blockSize + ii)*blockSize + jj]);
            printf("\n");
        }
    }
#endif
#elif defined(USE_MKL) && !defined(USE_MPI)
    sparse_status_t err;
    sparse_matrix_t coordMatrix;

    mat->blockSize = blockSize;
    mat->i = MALLOC(mat->nonzeros * sizeof(*(mat->i)));
    mat->j = MALLOC(mat->nonzeros * sizeof(*(mat->j)));
    mat->value = MALLOC(mat->nonzeros * sizeof(*(mat->value)));
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
    FREE(mat->i);
    FREE(mat->j);
    FREE(mat->value);
    if((err = mkl_sparse_destroy(coordMatrix)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "%i:  failed %i\n", wrank, err);
        exit(60);
    }
#else
    assert(blockSize > 0); // to suppress compiler warning
    mat->i = MALLOC(mat->nonzeros * sizeof(*mat->i));
    mat->j = MALLOC(mat->nonzeros * sizeof(*mat->j));
    mat->value = MALLOC(mat->nonzeros * sizeof(*mat->value));
    for(k=0; k<mat->nonzeros; k++){
        nread = fscanf(fp, "%u%u%le", mat->i+k, mat->j+k,
                       mat->value+k);
        mat->i[k] -= 1; mat->j[k] -= 1; // switch to zero-based indexing
        if(nread < 3) {
            fprintf(stderr, "%i:  Error reading %s, element %i\n",
                    wrank, fileName, k);
            break;
        }
    }
#endif

    /* As a first step, list range of matrix elements
       local to this process.  */
#ifdef USE_MPI
    mat->lowerRow = blockSize*rankIndex(wrank, wsize, mat->rows/blockSize);
    mat->lowerColumn = blockSize*
        rankIndex(wrank, wsize, mat->columns/blockSize);
#define MAX(A, B) (A>B?A:B)
    mat->gather = malloc(MAX(mat->rows, mat->columns)*sizeof(*mat->gather));
#endif
}

void sparseMatrixFree(SparseMatrix *mat) {
#if defined(USE_MKL) && !defined(USE_BLOCK) && !defined(USE_MPI)
    mkl_sparse_destroy(mat->a);
#else
    FREE(mat->value);
    FREE(mat->i);
    FREE(mat->j);
#endif
#ifdef USE_MPI
    FREE(mat->gather);
#endif
}

/* in and out must be distinct */
void matrixVector(const SparseMatrix *a,
                  const mat_int lin, const doublereal *in,
                  const mat_int lout, doublereal *out) {
#ifdef USE_MPI
    int wrank, wsize;
    MPI_Comm mpicom = a->mpicom;
    MPI_Comm_rank(mpicom, &wrank);
    MPI_Comm_size(mpicom, &wsize);
#else
    assert(lin == a->columns);
    assert(lout == a->rows);
#endif
#ifdef USE_BLOCK
    mat_int i, j, k, lower, upper;
    double *matp = a->value;
    const integer n = a->blockSize, n2=n*n;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T', normal='N';
#ifdef USE_MPI
    int rank;
    double *gather = a->gather;
    mat_int offset = 0;

    memcpy(gather + a->lowerColumn, in, lin*sizeof(*gather));
    for(rank=0; rank<wsize; rank++) {
        j = lin;
        MPI_Bcast(&j, 1, _MPI_MAT_INT, rank, mpicom);
        if(rank==wrank)
            assert(offset == a->lowerColumn);
        MPI_Bcast(gather+offset, j, MPI_DOUBLE, rank, mpicom);
        offset = offset + j;
    }
    lower = a->lowerRow;
#else
    const double *gather = in;

    lower = 0;
#endif
    upper = lower + lout;

    memset(out, 0, lout * sizeof(double));
    for(k=0; k<a->blocks; k++, matp+=n2) {
        i = a->i[k];
        j = a->j[k];
        if(i >= lower && i < upper)
        DGEMV(&trans, &n, &n, &one,
              matp, &n, gather + j, &inc, &one,
              out + i - lower, &inc);

        if(a->descr == 's' && i != j && j >= lower && j < upper) {
            DGEMV(&normal, &n, &n, &one,
                  matp, &n, gather + i, &inc, &one,
                  out + j - lower, &inc);
        }
    }
#elif defined(USE_MKL) && !defined(USE_MPI)
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
    mat_int i, j, k, lower, upper;
    double value;

#ifdef USE_MPI
    int rank;
    double *gather = a->gather;
    mat_int offset = 0;

    memcpy(gather + a->lowerColumn, in, lin*sizeof(*gather));
    for(rank=0; rank<wsize; rank++) {
        j = lin;
        MPI_Bcast(&j, 1, _MPI_MAT_INT, rank, mpicom);
        if(rank==wrank)
            assert(offset == a->lowerColumn);
        MPI_Bcast(gather+offset, j, MPI_DOUBLE, rank, mpicom);
        offset = offset + j;
    }
    lower = a->lowerRow;
#else
    const double *gather = in;

    lower = 0;
#endif
    upper = lower + lout;

    memset(out, 0, lout * sizeof(*out));
    for(k=0; k<a->nonzeros; k++) {
        i = a->i[k];
        j = a->j[k];
        value = a->value[k];
        if(i >= lower && i < upper)
            out[i-lower] += value * gather[j];
        if(a->descr == 's' && j != i && j >= lower && j < upper) {
                out[j-lower] += value * gather[i];
        }
    }
#endif
}

void vectorMatrix(const SparseMatrix *a,
                  const mat_int lin, const doublereal *in,
                  const mat_int lout, doublereal *out) {
#ifdef USE_MPI
    int wrank, wsize;
    MPI_Comm mpicom = a->mpicom;
    MPI_Comm_rank(mpicom, &wrank);
    MPI_Comm_size(mpicom, &wsize);
#else
    assert(lin == a->rows);
    assert(lout == a->columns);
#endif
#ifdef USE_BLOCK
    mat_int i, j, k, lower, upper;
    double *matp = a->value;
    const integer n = a->blockSize, n2 = n*n;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T', normal='N';

#ifdef USE_MPI
    int rank;
    double *gather = a->gather;
    mat_int offset = 0;

    memcpy(gather + a->lowerRow, in, lin*sizeof(*gather));
    for(rank=0; rank<wsize; rank++) {
        j = lin;
        MPI_Bcast(&j, 1, _MPI_MAT_INT, rank, mpicom);
        if(rank==wrank)
            assert(offset == a->lowerRow);
        MPI_Bcast(gather+offset, j, MPI_DOUBLE, rank, mpicom);
        offset = offset + j;
    }
    lower = a->lowerColumn;
#else
    const double *gather = in;

    lower = 0;
#endif
    upper = lower + lout;

    memset(out, 0, lout * sizeof(*out));
    for(k=0; k<a->blocks; k++, matp+=n2) {
        i = a->i[k];
        j = a->j[k];
        if(j >= lower && j < upper)
            DGEMV(&normal, &n, &n, &one,
                  matp, &n, gather + i, &inc, &one,
                  out + j - lower, &inc);

        if(a->descr == 's' && i != j && i >= lower && i < upper) {
            DGEMV(&trans, &n, &n, &one,
                  matp, &n, gather + j, &inc, &one,
                  out + i - lower, &inc);
        }
    }
#elif defined(USE_MKL) && !defined(USE_MPI)
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
    mat_int i, j, k, lower, upper;
    double value;
#ifdef USE_MPI
    int rank;
    mat_int offset = 0;
    double *gather = a->gather;

    memcpy(gather + a->lowerRow, in, lin*sizeof(*gather));
    for(rank=0; rank<wsize; rank++) {
        j = lin;
        MPI_Bcast(&j, 1, _MPI_MAT_INT, rank, mpicom);
        if(rank==wrank)
            assert(offset == a->lowerRow);
        MPI_Bcast(gather+offset, j, MPI_DOUBLE, rank, mpicom);
        offset = offset + j;
    }
    lower = a->lowerColumn;
#else
    const double *gather = in;

    lower = 0;
#endif
    upper = lower + lout;

    memset(out, 0, lout * sizeof(*out));
    for(k=0; k<a->nonzeros; k++) {
        i = a->i[k];
        j = a->j[k];
        value = a->value[k];
        if(j >= lower && j < upper)
            out[j-lower] += value * gather[i];
        if(a->descr == 's' && i != j && i >= lower && i < upper) {
            out[i-lower] += value * gather[j];
        }
    }
#endif
}
