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
#include "mmio.h"


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


void sparseMatrixRead(SparseMatrix *mat, char *fileName, char descr,
                      int tFlag, int blockSize, int chunkSize,
                      _MPI_Comm mpicom) {
    mat_int i, j, k;
#if defined(USE_BLOCK) || defined(USE_MPI) || !defined(USE_MKL)
    mat_int l, lowerCount, upperCount, nz;
#endif
    FILE *fp;
    MM_typecode matcode;
    int nread, wrank;
#ifdef USE_MPI
    mat_int lower, upper;
    int wsize;

    mat->mpicom = mpicom;
    MPI_Comm_rank(mpicom, &wrank);
    MPI_Comm_size(mpicom, &wsize);
#else
    wrank = 0;
    assert(mpicom == NULL);
#endif


    /*
           Read Matrix Market file into SparseMatrix
    */
    if(wrank ==0)
        printf("Opening file %s", fileName);
    if((fp = fopen(fileName, "r")) == NULL) {
        fprintf(stderr, "%i:  Could not open file %s\n",
                wrank, fileName);
        exit(44);
    }
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr, "%i:  Could not process Matrix Market banner.\n",
                wrank);
        exit(701);
    }
    if(descr == 'g' && !mm_is_general(matcode)) {
        fprintf(stderr, "%i:  Gauge matrix should general\n", wrank);
        fprintf(stderr, "%i:  Market Market type: [%s]\n",
                wrank, mm_typecode_to_str(matcode));
        exit(702);
    }
    if(!(mm_is_real(matcode) && mm_is_matrix(matcode) && 
         mm_is_coordinate(matcode))) {
        fprintf(stderr, "%i:  Matrix should be real, coordinate.\n",
                wrank);
        fprintf(stderr, "%i:  Market Market type: [%s]\n",
                wrank, mm_typecode_to_str(matcode));
        exit(702);
    }
    if (mm_read_mtx_crd_size(fp,
                             tFlag?&mat->columns:&mat->rows,
                             tFlag?&mat->rows:&mat->columns,
                             &mat->nonzeros) !=0) {
        fprintf(stderr, "%i:  Cannot read dimensions\n", wrank);
        exit(703);
    }
    if(wrank ==0)
        printf(" with %i nonzero elements.\n", mat->nonzeros);
    mat->i = MALLOC(mat->nonzeros * sizeof(*mat->i));
    mat->j = MALLOC(mat->nonzeros * sizeof(*mat->j));
    mat->value = MALLOC(mat->nonzeros * sizeof(*mat->value));
    for(k=0; k<mat->nonzeros; k++){
        nread = fscanf(fp, "%u%u%le", tFlag?&j:&i,
                       tFlag?&i:&j, mat->value+k);
         // switch to zero-based indexing, possible type conversion
        mat->i[k] = i - 1; mat->j[k] = j - 1;
        if(nread < 3) {
            fprintf(stderr, "%i:  Error reading %s, element %i\n",
                    wrank, fileName, k);
            break;
        }
    }
    fclose(fp);


    /*
       Define quantities needed by matrixVector
    */
#ifdef USE_MPI
    mat->lowerColumn = blockSize*
        rankIndex(wrank, wsize, mat->columns/blockSize);
#define MAX(A, B) (A>B?A:B)
    mat->gather = malloc(MAX(mat->rows, mat->columns)*
                         sizeof(*mat->gather));
    mat->localTime = 0.0;
    mat->mpiTime = 0.0;
    mat->count = 0;
#endif


    /*
       Fill the other triangle of a symmetric matrix

       This assumes the matrix has no block structure assigned
    */
#if defined(USE_BLOCK) || defined(USE_MPI) || !defined(USE_MKL)
    if(descr == 's') {
        lowerCount = 0; upperCount = 0;
        for(k=0; k<mat->nonzeros; k++) {
            if(mat->i[k] > mat->j[k])
                lowerCount++;
            if(mat->i[k] < mat->j[k])
                upperCount++;
        }
        // Verify only one triangle is filled
        assert(lowerCount == 0 || upperCount == 0);
        nz = mat->nonzeros + lowerCount + upperCount;
        mat->i = REALLOC(mat->i, nz * sizeof(*mat->i));
        mat->j = REALLOC(mat->j, nz * sizeof(*mat->j));
        mat->value = REALLOC(mat->value, nz * sizeof(*mat->value));
        for(k=0, l=mat->nonzeros; k<mat->nonzeros; k++)
            if(mat->i[k] != mat->j[k]) {
                mat->i[l] = mat->j[k];
                mat->j[l] = mat->i[k];
                mat->value[l] = mat->value[k];
                l++;
            }
        assert(l == nz);
        mat->nonzeros = nz;
    }
#endif


    /*
      Sort elements into ascending rows

      This is needed by the MKL matrices as well
      as the block creation algorithm.

      Row order is lost in the case of a transpose or 
      filling of a symmetric matrix
    */
    if(tFlag || descr=='s') {
#ifdef USE_BLOCK
        // Current matrix structure (needed for sort routine)
        mat->blockSize = 1;
        mat->blocks = mat->nonzeros;
#endif
        sortMatrix(mat, 1);
#if 1  // Test that the rows are in order.
        for(k=0, i=0, j=0; k<mat->nonzeros && j<10; k++) {
            if((mat_int) mat->i[k]<i) {
                if(j==0)
                    fprintf(stderr, "%i:  Sorted matrix %s, "
                            "tFlag=%i, exceptions:\n",
                            wrank, fileName, tFlag);
                fprintf(stderr, "%i:    k=%u i=%u j=%u value=%e lasti=%i \n",
                        wrank, k, mat->i[k], mat->j[k], mat->value[k], i);
                j++;
            }
            i = mat->i[k];
        }
        if(j>0)
            exit(88);
#endif
    }


    /*
        In the MPI case, filter out rows belonging to this
        process and shift row indices.
    */
#ifdef USE_MPI
    lower = blockSize*rankIndex(wrank, wsize, mat->rows/blockSize);
    upper = lower + blockSize*localSize(wrank, wsize, mat->rows/blockSize);
    for(k=0, l=0; k<mat->nonzeros; k++)
        if(mat->i[k] >= lower && mat->i[k] < upper) {
            mat->i[l] = mat->i[k] - lower;
            mat->j[l] = mat->j[k];
            mat->value[l] = mat->value[k];
            l++;
        }
    mat->i = REALLOC(mat->i, l*sizeof(*mat->i));
    mat->j = REALLOC(mat->j, l*sizeof(*mat->j));
    mat->value = REALLOC(mat->value, l*sizeof(*mat->value));
    // Sanity test:  make sure no elements were dropped.
    MPI_Allreduce(&l, &nz, 1, _MPI_MAT_INT,
                  MPI_SUM, mpicom); 
    assert(nz == mat->nonzeros);
    mat->nonzeros = l;
#endif


#ifdef USE_BLOCK
    mat_int ii = -1, jj = -1, lastii = 0, blockCols = 0;
    mat_int *ip = mat->i, *jp = mat->j;
    mat_int *blockColp;
    const mat_int maxBlockCol = mat->columns/blockSize;
    double value = 0.0, *blockRowValue, *blockValuep;
    double *valuep = mat->value;
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
            ii = ip[j];
            jj = jp[j];
            value = valuep[j];
            if(ii<lastii)
                printf("%i:  matrix %s bad order ii=%i lastii=%i jj=%i j=%i nonzeros=%i\n",
                       wrank, fileName, ii, lastii, jj, j, mat->nonzeros);
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
        blockRowValue[(k*blockSize + ii%blockSize)*blockSize +
                      jj%blockSize] = value;
    }
    free(blockRowValue);
    free(blockColFlag);
    free(blockColp);
    FREE(ip);
    FREE(jp);
    FREE(valuep);
#if 0
    printf("Print matrix for %s, tFlag=%i:\n", fileName, tFlag);
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
    if(mat->blocks*mat->blockSize*mat->blockSize != mat->nonzeros)
        printf("block elements=%u nonzeros=%i\n",
               mat->blocks*mat->blockSize*mat->blockSize,
               mat->nonzeros);
#endif

    sortMatrix(mat, chunkSize);
#elif defined(USE_MKL) && !defined(USE_MPI)
    sparse_status_t err;
    sparse_matrix_t coordMatrix;

    mat->blockSize = blockSize;
    if(descr == 's') {
        mat->descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        mat->descr.mode = SPARSE_FILL_MODE_LOWER;
        mat->descr.diag = SPARSE_DIAG_NON_UNIT;
    } else if(descr == 'g')
        mat->descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    else {
        fprintf(stderr, "Unknown descr\n");
        exit(890);
    }
    /* Create MKL data structure */
    if((err=mkl_sparse_d_create_coo(&coordMatrix, SPARSE_INDEX_BASE_ZERO,
              mat->rows, mat->columns, mat->nonzeros,
                               mat->i, mat->j, mat->value)) !=
           SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_create_coo failed %i\n", err);
        exit(55);
    }
    if((err = mkl_sparse_convert_bsr(coordMatrix, blockSize,
                            SPARSE_LAYOUT_ROW_MAJOR,
                            SPARSE_OPERATION_NON_TRANSPOSE,
                            &mat->a)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_convert_bsr failed %i\n", err);
        exit(56);
    }
#if 0  // Returns unsupported error
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
        fprintf(stderr, "mkl_sparse_optimize failed %i\n", err);
        exit(59);
    }
    FREE(mat->i);
    FREE(mat->j);
    FREE(mat->value);
    if((err = mkl_sparse_destroy(coordMatrix)) != SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "destroy failed %i\n", err);
        exit(60);
    }
    assert(chunkSize > 0);  // Suppress compiler message
#else
    assert(blockSize > 0);  // Suppress compiler message
    sortMatrix(mat, chunkSize);
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
void matrixVector(SparseMatrix *a,
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
    mat_int i, j, k;
    double *matp = a->value;
    const integer n = a->blockSize, n2=n*n;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T';
#ifdef USE_MPI
    int rank;
    double *gather = a->gather, t0, t1;
    mat_int offset = 0;

    t0 = MPI_Wtime();
    memcpy(gather + a->lowerColumn, in, lin*sizeof(*gather));
    for(rank=0; rank<wsize; rank++) {
        j = lin;
        MPI_Bcast(&j, 1, _MPI_MAT_INT, rank, mpicom);
        if(rank==wrank)
            assert(offset == a->lowerColumn);
        MPI_Bcast(gather+offset, j, MPI_DOUBLE, rank, mpicom);
        offset = offset + j;
    }
    assert(offset == a->columns);
    a->mpiTime += (t1 = MPI_Wtime()) - t0;
#else
    const double *gather = in;
#endif

    memset(out, 0, lout * sizeof(double));
    for(k=0; k<a->blocks; k++, matp+=n2) {
        i = a->i[k];
        j = a->j[k];
        DGEMV(&trans, &n, &n, &one,
              matp, &n, gather + j, &inc, &one,
              out + i, &inc);
    }
#ifdef USE_MPI
    a->localTime += MPI_Wtime() - t1;
    a->count += 1;
#endif
#elif defined(USE_MKL) && !defined(USE_MPI)
    sparse_status_t err;
    const double alpha = 1.0, beta=0.0;
    double d;

    if((err=mkl_sparse_d_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
                  alpha, a->a, a->descr, in, beta, out, &d)) !=
                              SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_dotmv failed %i\n", err);
#include <signal.h>
        raise(SIGINT);
        exit(69);
    }
#else
    mat_int i, j, k;
    double value;

#ifdef USE_MPI
    int rank;
    double *gather = a->gather, t0, t1;
    mat_int offset = 0;

    t0 = MPI_Wtime();
    memcpy(gather + a->lowerColumn, in, lin*sizeof(*gather));
    for(rank=0; rank<wsize; rank++) {
        j = lin;
        MPI_Bcast(&j, 1, _MPI_MAT_INT, rank, mpicom);
        if(rank==wrank) {
            assert(offset == a->lowerColumn);
        }
        MPI_Bcast(gather+offset, j, MPI_DOUBLE, rank, mpicom);
        offset = offset + j;
    }
    a->mpiTime += (t1 = MPI_Wtime()) - t0;
#else
    const double *gather = in;
#endif

    memset(out, 0, lout * sizeof(*out));
    for(k=0; k<a->nonzeros; k++) {
        i = a->i[k];
        j = a->j[k];
        value = a->value[k];
        out[i] += value * gather[j];
    }
#ifdef USE_MPI
    a->localTime += MPI_Wtime() - t1;
    a->count += 1;
#endif
#endif
}
