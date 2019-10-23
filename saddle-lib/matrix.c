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


void sparseMatrixRead(SparseMatrix *mat, char *fileName, char descr,
                      int tFlag, int blockSize, int chunkSize,
                      _MPI_Comm mpicom) {
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
    int p, q, r, jrank, lastJrank;
    mat_int j0;
    TaskList *task;
#ifdef USE_BLOCK
    const mat_int b1 = mat->blockSize, blocks = mat->blocks;
#else
    const mat_int b1 = 1, blocks = mat->nonzeros;
#endif

    for(p=0; p<2; p++)
        mat->gather[p] = MALLOC(b1*maxLocalSize(wsize, mat->columns/b1)*
                                sizeof(**mat->gather));

    sortMatrixLocal(mat, wrank, wsize);

    /* verify order, mark beginning and end of
       each process block, and shift column indices.
       Create row pointers. */
    mat->task = malloc(wsize*sizeof(TaskList));
    mat->task[0].start = 0;
    for(k=0, mat->taskCount=0; ; k++) {
        if(k < blocks) {
            jrank = indexRank(mat->j[k]/b1, wsize, mat->columns/b1);
            assert(k==0 || jrank >= lastJrank);
            if(k==0 || jrank != lastJrank)
                j0 = b1*rankIndex(jrank, wsize, mat->columns/b1);
            mat->j[k] -= j0;
        }
        if(k==blocks || (k>0 && jrank > lastJrank)) {
            // End of process block
            mat->task[mat->taskCount].doRank = lastJrank;
            mat->task[mat->taskCount].doSize =
                b1*localSize(lastJrank, wsize, mat->columns/b1);
            mat->task[mat->taskCount].end = k;
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
        if(p < mat->taskCount)
            needs=mat->task[p].doRank;
        else
            needs=-1;
        mat->task[p].sendTo = -1; // Don't send is default.
        // Start receive for next round of work
        if(p+1 < mat->taskCount && mat->task[p+1].doRank != wrank) {
            mat->task[p].receiveFrom = mat->task[p+1].doRank;
            mat->task[p].receiveSize = mat->task[p+1].doSize;
        } else {
            mat->task[p].receiveFrom = -1;
            mat->task[p].receiveSize = 0;
        }
#ifdef USE_MPI
        MPI_Allgather(&needs, 1, MPI_INT, allNeeds, 1, MPI_INT, mpicom);
#else
        allNeeds[0] = needs;
#endif
        for(q = 0; q<wsize; q++) {
            if(allNeeds[q] == wrank && q!= wrank) {
                /* True if each row starts with the
                   process-local block and the process-local
                   block is non-empty. */
                assert(p>0);
                /* Demand that the work is evenly distributed 
                   across processes.  See comment below. */
                assert(p-1 < mat->taskCount);
                /* Only one send per task. True if the process-block
                 structure is the same for each process. If needed, 
                 one could extend this to the general case.

                But, in that case, the work on each process must be 
                fairly asymmetric and our overall strategy may not
                be efficient anyway. */
                assert(mat->task[p-1].sendTo == -1);
                mat->task[p-1].sendTo = i;
            }
        }
    }
    free(allNeeds);

    /* print out the results */
    for(p = 0; p<wsize; p++) {
        if(wrank == 0)
            printf("task columns:  start, end, doRank, doSize, sendTo, "
                  "receiveFrom, receiveSize\n");
#ifdef USE_MPI
        MPI_Barrier(mpicom);
#endif
        if(p == wrank) {
            printf("%i:  Task list, taskCount=%i\n",
                   wrank, mat->taskCount);
            for(q=0; q<wsize; q++) {
                task = mat->task + q;
                printf("%i:  %i %i %i %i %i %i %i\n", wrank,
                    task->start, task->end, task->doRank, task->doSize,
                    task->sendTo, task->receiveFrom, task->receiveSize);
            }
        }
    }
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
    int p, wsize;
    mat_int i, j, k, start, end;
#ifdef USE_MPI
    double *gather;
    int wrank;
    MPI_Comm mpicom = a->mpicom;
    MPI_Comm_rank(mpicom, &wrank);
    MPI_Comm_size(mpicom, &wsize);
    double t0, t1;
#ifdef USE_TASK
    MPI_Request sent, received[2];
#else
    mat_int offset = 0;
    int rank;
#endif
#else // USE_MPI
    assert(lin == a->columns);
    assert(lout == a->rows);
    wsize = 1;
    const double *gather;
#endif
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


    memset(out, 0, lout * sizeof(double));
    for(p = 0; p < taskCount; p++) {
#ifdef USE_TASK
        task = a->task + p;
#endif

        /* Gather and send input vector */

#ifdef USE_MPI
        t0 = MPI_Wtime();
#ifdef USE_TASK
        if(task->sendTo != -1) {
            MPI_Isend(in, lin, MPI_DOUBLE,
               task->sendTo, MPI_ANY_TAG, mpicom, &sent);
        }
        if(task->receiveFrom != -1) {
            // Cycle through the two gather arrays.
            MPI_Irecv(a->gather[(p+1)%2], task->receiveSize, MPI_DOUBLE,
                      task->receiveFrom, MPI_ANY_TAG, mpicom,
                      received + ((p+1)%2));
        }
        if(task->doRank != wrank && task->doRank != -1)
            MPI_Wait(received+(p%2), MPI_STATUS_IGNORE);
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
        if(task->doRank)
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
        for(k=start, ap=a->value + start; k<end; k++, ap += b2) {
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
        a->localTime += MPI_Wtime() - t1;
#endif
        } // loop over tasks

#ifdef USE_MPI
    a->count += 1;
#endif
#endif  // MKL or not
}
