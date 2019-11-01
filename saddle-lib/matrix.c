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
// Lowest rank has maximum size
mat_int maxLocalSize(const int wsize, const mat_int n) {
    return localSize(0, wsize, n);
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


/*
   Debug print of matrix-vector product
*/
void testMatrixVector(SparseMatrix *mat, double *in) {
    double *y;
    mat_int k, nrow, ncol;
    int wrank, wsize, i;
#ifdef USE_MPI
    MPI_Comm_size(mat->mpicom, &wsize);
    MPI_Comm_rank(mat->mpicom, &wrank);
#else
    wsize = 1;
    wrank = 0;
#endif
#ifdef USE_BLOCK
    const int b1 = mat->blockSize;
#else
    const int b1 = 1;
#endif

    nrow = b1*localSize(wrank, wsize, mat->rows/b1);
    ncol = b1*localSize(wrank, wsize, mat->columns/b1);
    y = malloc(nrow * sizeof(double));
    matrixVector(mat, ncol, in, nrow, y);
    for(i=0; i<wsize; i++) {
        if(i == wrank) {
            if(wrank == 0)
                printf("Matrix-vector product result:\n");
            for(k = 0; k<nrow; k++) {
            printf("%i:  %le (%i)\n", wrank, y[k], k);
            }
        }
#ifdef USE_MPI
        MPI_Barrier(mat->mpicom);
#endif
    }
    free(y);
}

void printMatrixBlocks(SparseMatrix *mat, int wrank, int wsize,
                       mat_int blockSize, char *fileName, int tFlag,
                       _MPI_Comm mpicom) {
    int rank;
    mat_int k, ii, jj;
#ifdef USE_BLOCK
    mat_int n = mat->blocks;
#else
    mat_int n = mat->nonzeros;
#endif
    for(rank=0; rank<wsize; rank++) {
#ifdef USE_MPI
        MPI_Barrier(mpicom);
#else
	assert(mpicom == NULL);
#endif
        if(rank==wrank) {
            printf("%i:  Print matrix for %s, tFlag=%i:\n",
                   wrank, fileName, tFlag);
            for(k=0; k<n; k++) {
                printf("%i:  %i %i  ", wrank, mat->i[k], mat->j[k]);
                for(ii=0; ii<blockSize; ii++) {
                    if(ii>0)
                        printf("%i:      ", wrank);
                    for(jj=0; jj<blockSize; jj++)
                        printf(" %.3f", mat->value[(k*blockSize + ii)*blockSize + jj]);
                    printf("\n");
                }
            }
#ifdef USE_BLOCK
            if(mat->blocks*mat->blockSize*mat->blockSize != mat->nonzeros)
                printf("block elements=%u nonzeros=%i\n",
                       mat->blocks*mat->blockSize*mat->blockSize,
                       mat->nonzeros);
#endif
        }
    }
}

void printTasks(SparseMatrix *mat, int wsize, int wrank, _MPI_Comm mpicom) {
#ifdef USE_TASK
    int p, q;
    TaskList *task;
    for(p = 0; p<wsize; p++) {
#ifdef USE_MPI
        MPI_Barrier(mpicom);
#else
	assert(mpicom == NULL);
#endif
        if(p == wrank) {
            if(wrank == 0)
                printf("task columns:\n    task, start, end, "
                       "doRank, sendTo, receiveFrom, receiveSize\n");
            for(q=0; q<mat->taskCount; q++) {
                task = mat->task + q;
                printf("%i:  %i %i %i %i %i %i %i\n", wrank, q,
                    task->start, task->end, task->doRank,
                       task->sendTo, task->receiveFrom, task->receiveSize);
            }
        }
    }
#endif
}


void sparseMatrixRead(SparseMatrix *mat, char *fileName, char descr,
                      int tFlag, int blockSize, int chunkSize,
                      int debug, _MPI_Comm mpicom) {
    mat_int i, j, k;
#if !USE_MKL_MATRIX
    mat_int l, lowerCount, upperCount, nz;
#endif
    FILE *fp;
    MM_typecode matcode;
    int nread, wrank, wsize;
#ifdef USE_MPI
    mat_int lower, upper;

    mat->mpicom = mpicom;
    MPI_Comm_rank(mpicom, &wrank);
    MPI_Comm_size(mpicom, &wsize);
#else
    wrank = 0;
    wsize = 1;
    assert(mpicom == NULL);
#endif

    // Currently this quantity is unused.
    assert(chunkSize == 1);

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
       Fill the other triangle of a symmetric matrix

       This assumes the matrix has no block structure assigned
    */
#if !USE_MKL_MATRIX
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

      Row order has been lost in the case of a transpose or 
      filling of a symmetric matrix

      Row (at least) ordering is needed by the MKL matrices as well
      as the block creation algorithm.
    */
    if(tFlag || descr=='s') {
#ifdef USE_BLOCK
        // Current matrix structure (needed for sort routine)
        mat->blockSize = 1;
        mat->blocks = mat->nonzeros;
#endif
        sortMatrixChunks(mat, 1);
#if 1  // Verify that the rows are in order.
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

        In the non-MPI case, this should do nothing.
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
    // Sanity test:  verify that no elements were dropped.
    MPI_Allreduce(&l, &nz, 1, _MPI_MAT_INT,
                  MPI_SUM, mpicom); 
    assert(nz == mat->nonzeros);
    if(debug>0) {
        printf("%i:  local rows from %i to %i\n", wrank,
            mat->rows, upper-lower);
        printf("%i:  nonzeros from %i to %i\n", wrank, mat->nonzeros, l);
    }
    mat->nonzeros = l;
#endif


    /* 
       Divide matrix into blocks or create MKL block matrix 
    */
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
    if(debug>0)
        printf("%i:  Create block matrix, blocks=%i, blockSize=%i from %i nonzeros\n",
            wrank, mat->blocks, mat->blockSize, mat->nonzeros);
    if(debug>1)
        printMatrixBlocks(mat, wrank, wsize, blockSize,
                          fileName, tFlag, mpicom);

#elif USE_MKL_MATRIX
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
        exit(54);
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
#else
    assert(blockSize > 0);  // Suppress compiler message
#endif


    /*
      Sort matrix into local process blocks
      Create schedule and shift column labels.

      Currently, as defined in Mathematica, the basis 
      is ordered according to one of the lattice dimensions.
      If the number of processes is less than the 
      length of the lattice, then the process of rank r will
      only have nonzero matrix elements with processes rank 
      (r-1)%wsize and (r+1)%wsize.

      In the non-MPI case, this should work trivially.
      In particular, the sort should have no effect.
    */
#ifdef USE_TASK
    int needs, *allNeeds = malloc(wsize*sizeof(*allNeeds));
    int p, q, jrank = -1, lastJrank = -1;
    mat_int j0 = 0;
#ifdef USE_BLOCK
    const mat_int b1 = mat->blockSize, blocks = mat->blocks;
#else
    const mat_int b1 = 1, blocks = mat->nonzeros;
#endif

    for(p=0; p<2; p++)
        mat->gather[p] = MALLOC(b1*maxLocalSize(wsize, mat->columns/b1)*
                                sizeof(*mat->value));

    /* For the debug print, processes take turns sorting the matrix.
       In the debug print case, change 2 flags in sort.c, examine
       results in Mathematica notebook test-sort.nb. */
    if(debug>1) { // debug print
        for(p=0; p<wsize; p++) {
#ifdef USE_MPI
            MPI_Barrier(mpicom);
#endif
            if(p==wrank) {
                printf("%i: Matrix %s%s sort:\n", wrank, fileName, tFlag?"T":"");
                sortMatrixLocal(mat, wrank, wsize, 1);
            }
        }
    } else {
        // Normal case
        sortMatrixLocal(mat, wrank, wsize, 0);
    }
    if(debug>1) {
        printf("%i:  Block matrix after localSort\n", wrank);
        printMatrixBlocks(mat, wrank, wsize, blockSize,
                          fileName, tFlag, mpicom);
    }

    /* Verify order, mark beginning and end of
       each process block, and shift column indices for each block. */
    mat->task = malloc(wsize*sizeof(TaskList));
    mat->task[0].start = 0;
    for(k=0, mat->taskCount=0; ; k++) {
        if(k < blocks) {
            jrank = indexRank(mat->j[k]/b1, wsize, mat->columns/b1);
            /* Verify (jrank-wrank)%wsize is non-decreasing
               This is not, strictly speaking, necessary, but it
               does verify that sortMatrixLocal() is working as
               intended. */
            /* The C standard for mod of negative numbers is screwy:
               add wsize to get normal behavior. */
            if(k>0 && (wsize+jrank-wrank)%wsize <
               (wsize+lastJrank-wrank)%wsize) {
                fprintf(stderr, "%i:  Expect non-decreasing rank.\n", wrank);
                fprintf(stderr, "%i:  matrix %s%s k=%i jrank=%i lastJrank=%i\n",
                        wrank, fileName, tFlag?"T":"", k, jrank, lastJrank);
                exit(55);
            }
            if(k==0 || jrank != lastJrank)
                j0 = b1*rankIndex(jrank, wsize, mat->columns/b1);
            mat->j[k] -= j0;
        }
        if(k==blocks || (k>0 && jrank != lastJrank)) {
            // End of process block
            mat->task[mat->taskCount].doRank = lastJrank;
            mat->task[mat->taskCount].end = k;
            mat->task[mat->taskCount].sendTo = -1; // Don't send is default.
            mat->taskCount++;
            if(k == blocks)
                break;
            assert(mat->taskCount<wsize);
            mat->task[mat->taskCount].start = k;
        }
        lastJrank = jrank;
    }
    /* Share the data needed by each process. */
    for(p=0; p<wsize; p++){
        // Start receive for next round of work
        if(p+1 < mat->taskCount && mat->task[p+1].doRank != wrank) {
            mat->task[p].receiveFrom = mat->task[p+1].doRank;
            mat->task[p].receiveSize =
                b1*localSize(mat->task[p+1].doRank, wsize, mat->columns/b1);
        } else {
            mat->task[p].receiveFrom = -1;
            mat->task[p].receiveSize = 0;
        }
        if(p < mat->taskCount)
            needs=mat->task[p].doRank;
        else
            needs=-1;
        // printf("%i:  p=%i needs=%i\n", wrank, p, needs);
#ifdef USE_MPI
        MPI_Allgather(&needs, 1, MPI_INT, allNeeds, 1, MPI_INT, mpicom);
#else
        allNeeds[0] = needs;
#endif
        for(q = 0; q<wsize; q++) {
            if(allNeeds[q] == wrank && q != wrank) {
                /* Verify that each row starts with the
                   process-local block and the process-local
                   block is non-empty. */
                if(p == 0) {
                    fprintf(stderr, "%i:  matrix %s%s q=%i allNeeds[q]=%i\n",
                            wrank, fileName, tFlag?"T":"", q, allNeeds[q]);
                    exit(77);
                }
                /* Demand that the work is evenly distributed 
                   across processes.  See comment below. */
                assert(p-1 < mat->taskCount);
                /* Only one send per task. True if the process-block
                 structure is the same for each process. If needed, 
                 one could extend this to the general case.

                But, in that case, the work on each process must be 
                fairly asymmetric and our overall strategy may not
                be efficient anyway. */
                if(mat->task[p-1].sendTo != -1) {
                    fprintf(stderr, "%i: matrix %s%s multiple sends.\n",
                            wrank, fileName, tFlag?"T":"");
                    fprintf(stderr, "%i: task %i sendTo=%i and %i.  "
                            "Maybe too many processes?\n",
                            wrank, p-1, mat->task[p-1].sendTo, q);
                    exit(79);
                }
                mat->task[p-1].sendTo = q;
            }
        }
    }
    free(allNeeds);

    if(debug>1)
        printMatrixBlocks(mat, wrank, wsize, blockSize,
                          fileName, tFlag, mpicom);
    if(debug>0)
        printTasks(mat, wsize, wrank, mpicom);
#else
    mat->gather = MALLOC(mat->columns*sizeof(*mat->gather));
    mat->lowerColumn = blockSize*
        rankIndex(wrank, wsize, mat->columns/blockSize);
#endif


    /* 
       Initialize profiling associated with MPI
    */
#ifdef USE_MPI
    mat->localTime = 0.0;
    mat->mpiTime = 0.0;
    mat->count = 0;
#endif
}

void sparseMatrixFree(SparseMatrix *mat) {
#if USE_MKL_MATRIX
    mkl_sparse_destroy(mat->a);
#else
    FREE(mat->value);
    FREE(mat->i);
    FREE(mat->j);
#endif
#ifdef USE_TASK
    int i;
    for(i=0; i<2; i++)
        FREE(mat->gather[i]);
    free(mat->task);
#else
    FREE(mat->gather);
#endif
}

/* in and out must be distinct */
void matrixVector(SparseMatrix *a,
                  const mat_int lin, const doublereal *in,
                  const mat_int lout, doublereal *out) {
#if USE_MKL_MATRIX
    sparse_status_t err;
    const double alpha = 1.0, beta=0.0;
    double d;

    assert(lin == a->columns);
    assert(lout == a->rows);
    if((err=mkl_sparse_d_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
                  alpha, a->a, a->descr, in, beta, out, &d)) !=
       SPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "mkl_sparse_d_dotmv failed %i\n", err);
        exit(69);
    }
#else // USE_MKL_MATRIX
    int p;
    mat_int i, j, k, start, end;
#ifdef USE_MPI
    int wrank;
    MPI_Comm mpicom = a->mpicom;
    MPI_Request sent = MPI_REQUEST_NULL, received = MPI_REQUEST_NULL;
    double t0, t1;
#ifdef USE_TASK
    const double *gather;
#else
    mat_int offset = 0;
    int rank, wsize;
    double *gather;
    MPI_Comm_size(mpicom, &wsize);
#endif
#else // USE_MPI
    const double *gather;
    assert(lin == a->columns);
    assert(lout == a->rows);
#endif // USE_MPI
#ifdef USE_TASK
    TaskList *task;
    int taskCount = a->taskCount;
#else
    int taskCount = 1;
#endif
#ifdef USE_BLOCK
    double *ap;
    const integer b1 = a->blockSize, b2=b1*b1;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T';
#endif


#ifdef USE_MPI
    MPI_Comm_rank(mpicom, &wrank);
    t0 = MPI_Wtime();
#endif
    memset(out, 0, lout * sizeof(double));
    for(p = 0; p < taskCount; p++) {
#ifdef USE_TASK
        task = a->task + p;
#endif

        /* Gather and send input vector */

#ifdef USE_MPI
#ifdef USE_TASK
        if(task->sendTo != -1) {
            MPI_Isend(in, lin, MPI_DOUBLE,
               task->sendTo, p+1, mpicom, &sent);
        }
        if(task->receiveFrom != -1) {
            // Alternate between the two gather arrays.
            MPI_Irecv(a->gather[(p+1)%2], task->receiveSize, MPI_DOUBLE,
                      task->receiveFrom, p+1, mpicom, &received);
        }
        if(task->doRank == wrank)
            gather = in;
        else if(task->doRank == -1)
            gather = NULL;
        else
            gather = a->gather[p%2];
#else // USE_TASK
        /* The strategy used here is very inefficient.
           It is intended only for debugging. */
        gather = a->gather;
        memcpy(gather + a->lowerColumn, in, lin*sizeof(*gather));
        for(rank=0; rank<wsize; rank++) {
            j = lin;
            MPI_Bcast(&j, 1, _MPI_MAT_INT, rank, mpicom);
            if(rank==wrank) {
#if 0
                printf("%i:  offset=%i lowerColumn=%i\n",
                    wrank, offset, a->lowerColumn);
#endif
                assert(offset == a->lowerColumn);
            }
            MPI_Bcast(gather+offset, j, MPI_DOUBLE, rank, mpicom);
            offset = offset + j;
        }
        assert(offset == a->columns);
#endif  // USE_TASK
        a->mpiTime += (t1=MPI_Wtime()) - t0;
#else // USE_MPI
        gather = in;
#endif // USE_MPI

        /* For a given task, do matrix-vector multiply */

#ifdef USE_TASK
        if(task->doRank == -1)
            continue;
        start = task->start;
        end = task->end;
#else // USE_TASK
        start = 0;
#ifdef USE_BLOCK
        end = a->blocks;
#else
        end = a->nonzeros;
#endif // USE_BLOCK
#endif // USE_TASK
#ifdef USE_BLOCK
        for(k=start, ap=a->value + b2*start; k<end; k++, ap += b2) {
            i = a->i[k];
            j = a->j[k];
            DGEMV(&trans, &b1, &b1, &one,
                  ap, &b1, gather + j, &inc, &one,
                  out + i, &inc);
        }
#else
        for(k=start; k<end; k++) {
            i = a->i[k];
            j = a->j[k];
            out[i] += a->value[k] * gather[j];
        }
#endif
#ifdef USE_MPI
        a->localTime += (t0 = MPI_Wtime()) - t1;

        /* Verify that any current messages are complete
           before starting the next task or exiting.
	   Subsequent code may modify the 'in' array. */
        MPI_Wait(&received, MPI_STATUS_IGNORE);
	MPI_Wait(&sent, MPI_STATUS_IGNORE);
#endif
        } // loop over tasks


#ifdef USE_MPI
    a->mpiTime += MPI_Wtime() - t0;
    a->count += 1;
#endif
#endif  // MKL or not
}
