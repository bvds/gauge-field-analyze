#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include "shifts.h"
#ifdef USE_MKL
#include "mkl_spblas.h"
#endif
#include "mmio.h"


/*
  Divide vectors and matrices up for mpi.  
  Without mpi, wrank = 0 and wsize = 1 gives
  the desired single-process behavior.

  "partitions" is used to specify how matrix rows
  are assigned to each mpi process.  For a matrix with n rows,
  each processor assigned an integer multiple of n/partitions.
  Since each proccess is also assigned some integral number 
  of blocks, n/partitions must be divisible by blockSize.
  
  The routine localsize(...) determines this multiple.
  The choice of "partitions" can be used to minimize the mpi 
  communication overhead.

  Setting partitions = n/blockSize will distribute blocks 
  evenly among the mpi processes.
 */
int indexRank(const mat_int i, const int wsize, const mat_int partitions) {
    if(i/(partitions/wsize + 1) < partitions%wsize)
        return i/(partitions/wsize + 1);
    else
        return (i-partitions%wsize)/(partitions/wsize);
}
mat_int localSize(const unsigned int wrank, const int wsize,
                  const mat_int partitions) {
    return partitions/wsize + ((wrank<(partitions%wsize))?1:0);
}
// Lowest rank has maximum size
mat_int maxLocalSize(const int wsize, const mat_int partitions) {
    return localSize(0, wsize, partitions);
}
mat_int rankIndex(const unsigned int wrank, const int wsize,
                  const mat_int partitions) {
    if(wrank < partitions%wsize)
        return wrank*(partitions/wsize + 1);
    else
        return wrank*(partitions/wsize) + partitions%wsize;
}


// Sanity tests for localSize, rankIndex, and indexRank
void rankSanityTest(mat_int partitions) {
    int i, wsize = 1;
    mat_int k;
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
#endif
    k = 0;
    for(i=0; i<wsize; i++) {
        assert(i == indexRank(rankIndex(i, wsize, partitions),
                              wsize, partitions));
        assert(localSize(i, wsize, partitions) ==
               rankIndex(i + 1, wsize, partitions)
               - rankIndex(i, wsize, partitions));
        k += localSize(i, wsize, partitions);
    }
    assert(k == partitions);
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

    nrow = (mat->rows/mat->partitions)*localSize(wrank, wsize, mat->partitions);
    ncol = (mat->columns/mat->partitions)*localSize(wrank, wsize, mat->partitions);
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
                       char *fileName, int tFlag) {
    int rank;
    mat_int k, ii, jj;
    mat_int n = mat->blocks, b1 = mat->blockSize;
    for(rank=0; rank<wsize; rank++) {
#ifdef USE_MPI
        MPI_Barrier(mat->mpicom);
#endif
        if(rank==wrank) {
            printf("%i:  Print matrix for %s, tFlag=%i:\n",
                   wrank, fileName, tFlag);
            for(k=0; k<n; k++) {
                printf("%i:  %i %i  ", wrank, mat->i[k], mat->j[k]);
                for(ii=0; ii<b1; ii++) {
                    if(ii>0)
                        printf("%i:      ", wrank);
                    for(jj=0; jj<b1; jj++)
                        printf(" %.3f", mat->value[(k*b1 + ii)*b1 + jj]);
                    printf("\n");
                }
            }
        }
    }
}

#ifdef USE_TASK
void printTasks(SparseMatrix *mat, int wsize, int wrank) {
    int p, q;
    TaskList *task;
    for(p = 0; p<wsize; p++) {
#ifdef USE_MPI
        MPI_Barrier(mat->mpicom);
#endif
        if(p == wrank) {
            if(wrank == 0)
                printf("task columns:\n    task, doRank, start, end, "
                       "sendTo, receiveFrom, receiveSize\n");
            for(q=0; q<mat->taskCount; q++) {
                task = mat->task + q;
                printf("%2i:  %2i %2i %3i %3i %3i %3i %3i\n", wrank, q,
                       task->doRank, task->start, task->end,
                       task->sendTo, task->receiveFrom, task->receiveSize);
            }
        }
    }
}
#endif


void sparseMatrixRead(SparseMatrix *mat, char *fileName, char descr,
                      int tFlag, int blockSize, mat_int partitions,
                      int chunkSize, int debug, _MPI_Comm mpicom) {
    mat_int i, j, k;
#if !USE_MKL_MATRIX
    mat_int l, lowerCount, upperCount, blocks;
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

    mat->partitions = partitions;
    // Currently this quantity is unused.
    assert(chunkSize == 1);
    
    /*
           Read Matrix Market file into SparseMatrix
    */
    if(wrank ==0 && debug>2)
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
    mat->blockSize = 1;  // The Matrix Market file has no block structure
    if (mm_read_mtx_crd_size(fp,
                             tFlag?&mat->columns:&mat->rows,
                             tFlag?&mat->rows:&mat->columns,
                             &mat->blocks) != 0) {
        fprintf(stderr, "%i:  Cannot read dimensions\n", wrank);
        exit(703);
    }
    if(wrank ==0 && debug>2)
        printf(" with %i nonzero elements.\n", mat->blocks);
    mat->i = malloc(mat->blocks * sizeof(*mat->i));
    mat->j = malloc(mat->blocks * sizeof(*mat->j));
    mat->value = malloc(mat->blocks * sizeof(*mat->value));
    for(k=0; k<mat->blocks; k++){
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
    */
#if !USE_MKL_MATRIX
    if(descr == 's') {
        assert(mat->blockSize == 1);  // assume no block structure
        lowerCount = 0; upperCount = 0;
        for(k=0; k<mat->blocks; k++) {
            if(mat->i[k] > mat->j[k])
                lowerCount++;
            if(mat->i[k] < mat->j[k])
                upperCount++;
        }
        // Verify only one triangle is filled
        assert(lowerCount == 0 || upperCount == 0);
        blocks = mat->blocks + lowerCount + upperCount;
        mat->i = realloc(mat->i, blocks * sizeof(*mat->i));
        mat->j = realloc(mat->j, blocks * sizeof(*mat->j));
        mat->value = realloc(mat->value, blocks * sizeof(*mat->value));
        for(k=0, l=mat->blocks; k<mat->blocks; k++)
            if(mat->i[k] != mat->j[k]) {
                mat->i[l] = mat->j[k];
                mat->j[l] = mat->i[k];
                mat->value[l] = mat->value[k];
                l++;
            }
        assert(l == blocks);
        if(debug>2 && wrank==0)
            printf("Fill symmetric matrix, blocks %i -> %i\n",
                   mat->blocks, blocks);
        mat->blocks = blocks;
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
        sortMatrixChunks(mat, 1);
#if 1  // Verify that the rows are in order.
        for(k=0, i=0, j=0; k<mat->blocks && j<10; k++) {
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
        process and shift row indices.  Even when blocks have
        not been assigned, make sure blocks are not divided 
        between processes.

        In the non-MPI case, this should do nothing.
    */
#ifdef USE_MPI
    mat_int rowPart = mat->rows/mat->partitions;
    lower = rowPart*rankIndex(wrank, wsize, mat->partitions);
    upper = lower + rowPart*localSize(wrank, wsize, mat->partitions);
    // One could expand this to handle the general case.
    assert(mat->blockSize == 1);
    for(k=0, l=0; k<mat->blocks; k++)
        if(mat->i[k] >= lower && mat->i[k] < upper) {
            mat->i[l] = mat->i[k] - lower;
            mat->j[l] = mat->j[k];
            mat->value[l] = mat->value[k];
            l++;
        }
    mat->i = realloc(mat->i, l*sizeof(*mat->i));
    mat->j = realloc(mat->j, l*sizeof(*mat->j));
    mat->value = realloc(mat->value, l*sizeof(*mat->value));
    // Sanity test:  verify that no elements were dropped.
    MPI_Allreduce(&l, &blocks, 1, _MPI_MAT_INT,
                  MPI_SUM, mpicom); 
    assert(blocks == mat->blocks);
    if(debug>2) {
        printf("%i:  local rows from %i to %i\n", wrank,
            mat->rows, upper-lower);
        printf("%i:  nonzeros from %i to %i\n", wrank, mat->blocks, l);
    }
    mat->blocks = l;
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
    assert(mat->blockSize == 1);
    assert(mat->rows%blockSize == 0);
    assert(mat->columns%blockSize == 0);

    mat->blockSize = blockSize;
    blocks = 0;
    blockRowValue = malloc(mat->columns*blockSize*sizeof(*blockRowValue));
    blockColFlag = malloc(maxBlockCol*sizeof(*blockColFlag));
    blockColp = malloc(maxBlockCol*sizeof(*blockColp));
    mat->i = malloc(0);
    mat->j = malloc(0);
    mat->value = malloc(0);
    memset(blockRowValue, 0, mat->columns*blockSize*sizeof(double));
    memset(blockColFlag, 0, maxBlockCol*sizeof(*blockColFlag));
    for(j=0;; j++){
        if(j<mat->blocks) {
            ii = ip[j];
            jj = jp[j];
            value = valuep[j];
            if(ii<lastii)
                printf("%i:  matrix %s bad order ii=%i lastii=%i jj=%i j=%i nonzeros=%i\n",
                       wrank, fileName, ii, lastii, jj, j, mat->blocks);
            assert(ii>=lastii);
        }
        if(j==mat->blocks || (lastii<ii && ii%blockSize == 0)) {
            // Retire block row and reset for a new one
            mat->i = realloc(mat->i,
                     (blocks + blockCols)*sizeof(*mat->i));
            mat->j = realloc(mat->j,
                     (blocks + blockCols)*sizeof(*mat->j));
            mat->value = realloc(mat->value,
                     (blocks + blockCols)*blockSize*blockSize*sizeof(double));
            for(i=0; i<blockCols; i++) {
                k = blockColp[i];
#if 0
                printf("  %i %i blocks=%i blockCols=%i\n",
                       j, k, blocks, blockCols);
#endif
                blockValuep = blockRowValue+k*blockSize*blockSize;
                memcpy(mat->value+blockSize*blockSize*blocks,
                       blockValuep, blockSize*blockSize*sizeof(double));
                memset(blockValuep, 0, blockSize*blockSize*sizeof(double));
                mat->i[blocks] = (lastii/blockSize)*blockSize;
                mat->j[blocks] = k*blockSize;
                blocks += 1;
                blockColFlag[k] = 0;
            }
            blockCols = 0;
        }
        if(j==mat->blocks)
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
    l = mat->blocks;
    mat->blocks = blocks;
    free(blockRowValue);
    free(blockColFlag);
    free(blockColp);
    free(ip);
    free(jp);
    free(valuep);
    if(debug>2)
        printf("%i:  Create block matrix, blocks=%i, blockSize=%i from %i nonzeros\n",
            wrank, mat->blocks, mat->blockSize, l);
    if(debug>3)
        printMatrixBlocks(mat, wrank, wsize, fileName, tFlag);

#elif USE_MKL_MATRIX
    sparse_status_t err;
    sparse_matrix_t coordMatrix;

    if(debug>2)
        printf("Create MKL BSR format matrix\n");
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
              mat->rows, mat->columns, mat->blocks,
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
    free(mat->i);
    free(mat->j);
    free(mat->value);
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
    const mat_int colPart = mat->columns/mat->partitions;
#ifdef USE_TASK
    int needs, *allNeeds = malloc(wsize*sizeof(*allNeeds));
    int p, q, jrank = -1, lastJrank = -1;
    mat_int j0 = 0;

    for(p=0; p<2; p++)
        MALLOC(mat->gather[p], colPart*maxLocalSize(wsize, mat->partitions)*
               sizeof(*mat->value));

    /* For the debug print, processes take turns sorting the matrix.
       In the debug print case, change 2 flags in sort.c, examine
       results in Mathematica notebook test-sort.nb. */
    if(debug>3) { // debug print
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
    if(debug>3) {
        printf("%i:  Block matrix after localSort\n", wrank);
        printMatrixBlocks(mat, wrank, wsize, fileName, tFlag);
    }

    /* Verify order, mark beginning and end of
       each process block, and shift column indices for each block.
       Also, create row pointers.

       iTemp handles cases where mat->iCount > k, which
       can happen, for instance, if the first blocks are zero.
    */
    mat_int *iTemp = malloc((mat->blocks + wsize)*sizeof(*iTemp));
    mat->ip = malloc((mat->blocks+1)*sizeof(*mat->ip));
    mat->task = malloc(wsize*sizeof(TaskList));
    mat->task[0].start = 0;
    for(k=0, mat->taskCount=0, mat->iCount=0; ; k++) {
        if(k < mat->blocks) {
            jrank = indexRank(mat->j[k]/colPart, wsize, mat->partitions);
            /* Verify (jrank-wrank)%wsize is non-decreasing.
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
                j0 = colPart*rankIndex(jrank, wsize, mat->partitions);
            mat->j[k] -= j0;  // shift j-index
        }
        if(k == mat->blocks || (k>0 && jrank != lastJrank)) {
            // End of process block
            mat->task[mat->taskCount].doRank = lastJrank;
            mat->task[mat->taskCount].end = mat->iCount;
            mat->task[mat->taskCount].sendTo = -1; // Don't send is default.
            mat->taskCount++;
            if(k == mat->blocks)
                break;
            assert(mat->taskCount < wsize);
            mat->task[mat->taskCount].start = mat->iCount;
        }
        // at this point, k<mat->blocks
        if(k==0 || mat->i[k] != mat->i[k-1] || jrank != lastJrank) {
            mat->ip[mat->iCount] = k;
            iTemp[mat->iCount] = mat->i[k];
            mat->iCount++;
        }
        lastJrank = jrank;
    }
    mat->i = realloc(mat->i, mat->iCount*sizeof(*mat->i));
    memcpy(mat->i, iTemp, mat->iCount*sizeof(*mat->i));
    free(iTemp);

    /* Share the data needed by each process. */
    for(p=0; p<wsize; p++){
        // Start receive for next round of work
        if(p+1 < mat->taskCount && mat->task[p+1].doRank != wrank) {
            mat->task[p].receiveFrom = mat->task[p+1].doRank;
            mat->task[p].receiveSize =
                colPart*localSize(mat->task[p+1].doRank, wsize, mat->partitions);
        } else {
            mat->task[p].receiveFrom = -1;
            mat->task[p].receiveSize = 0;
        }
        if(p < mat->taskCount)
            needs = mat->task[p].doRank;
        else
            needs = -1;
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
                    fprintf(stderr, "%i:  matrix %s%s, allNeeds[%i]=%i\n",
                            wrank, fileName, tFlag?"T":"", q, allNeeds[q]);
                    exit(77);
                }
                /* Demand that the work is fairly evenly distributed 
                   across processes.  See comment below. */
                assert(p-1 < mat->taskCount);
                /* Max of one send per task. True if the process-block
                 structure is fairly even across processes. If needed,
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

    if(debug>3)
        printMatrixBlocks(mat, wrank, wsize, fileName, tFlag);
    if(debug>2)
        printTasks(mat, wsize, wrank);
#else // USE_TASK
    MALLOC(mat->gather, mat->columns*sizeof(*mat->gather));
    mat->lowerColumn = colPart*rankIndex(wrank, wsize, mat->partitions);

    /* Create row pointers */
    mat->ip = malloc((mat->blocks+1)*sizeof(*mat->ip));
    for(k=0, mat->iCount=0; k<mat->blocks; k++) {
        if(k==0 || mat->i[k] != mat->i[k-1]) {
            assert(mat->iCount <= k); // Don't mangle mat->i
            mat->ip[mat->iCount] = k;
            mat->i[mat->iCount] = mat->i[k];
            mat->iCount++;
        }
    }
    mat->i = realloc(mat->i, mat->iCount*sizeof(*mat->i));
#endif // USE_TASK
    mat->ip[mat->iCount] = mat->blocks;
    mat->ip = realloc(mat->ip, (mat->iCount+1)*sizeof(*mat->ip));
    if(debug>3) {
        printf("%i:  Matrix blocks=%i\n", wrank, mat->blocks);
        for(k=0; k<mat->iCount; k++) {
            printf("%i:  k=%3i ip=%3i,%3i i=%3i\n",
               wrank, k, mat->ip[k], mat->ip[k+1], mat->i[k]);
        }
    }

    /* 
       Initialize profiling associated with MPI
    */
    mat->localTime = 0.0;
    mat->mpiTime = 0.0;
    mat->count = 0;
}


/*
  Print stats associated with matrix, so far
*/
void sparseMatrixStats(SparseMatrix *mat, char *matrixName) {
    int wrank;
    double nb, b1, opRate;
    double localTime, mpiTime;

#ifdef USE_MPI
    int wsize;
    MPI_Comm_size(mat->mpicom, &wsize);
    MPI_Comm_rank(mat->mpicom, &wrank);
#else
    wrank = 0;
#endif
    localTime = mat->localTime;
    mpiTime = mat->mpiTime;
    nb = mat->blocks;
    b1 = mat->blockSize;

#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &nb, 1, MPI_DOUBLE,
                  MPI_SUM, mat->mpicom);
    /*  Averaging is a bit silly
        Maybe finding the max and min would make more sense. */
    MPI_Allreduce(MPI_IN_PLACE, &localTime, 1, MPI_DOUBLE,
                  MPI_SUM, mat->mpicom);
    MPI_Allreduce(MPI_IN_PLACE, &mpiTime, 1, MPI_DOUBLE,
                  MPI_SUM, mat->mpicom);
    localTime /= wsize;
    mpiTime /= wsize;
#endif
    opRate = 1.0e-9*mat->count/localTime;
    if(wrank ==0 ) {
        printf("Matrix %s ops=%i, localTime=%.2f s, mpiTime=%.2f s\n",
               matrixName, mat->count, localTime, mpiTime);
               /*
                  For AMD EPYC 7301 16-Core Processor:
                  Theoretical max flopss is:
                     2.2GHz*16cores*8flops/core = 281 GFLOP/s
                  We count addition and multiplication as separate flops.

                  Theoretical memory bandwidth is 
                     2.666GHz*8bytes/channel*8channels = 170.62GB/s.
                  STREAM benchmark is memory reads plus writes.

                  According to the STREAM benchmark at
                  https://www.dell.com/support/article/us/en/04/sln313856/
                  we should expect roughly 126-130 GB/s for one socket
                  and 32 GB/s for one NUMA node.
                  They used 2400MHz memory, while I have 2666MHz memory.
               */
        printf("    global %.1f GFLOP/s, stream=%.1f GB/s\n",
               2*b1*b1*nb*opRate,  // multiply and add
	       /* optimal block compressed sparse row matrix storage
		  ignore row sums, since they are kept in cache */
	       opRate*(nb*(b1*(b1+1)*sizeof(double)+sizeof(mat_int)) +
                       mat->rows*sizeof(double)));
    }
}


void sparseMatrixFree(SparseMatrix *mat) {
#if USE_MKL_MATRIX
    mkl_sparse_destroy(mat->a);
#else
    free(mat->value);
    free(mat->i);
    free(mat->j);
    free(mat->ip);
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
    mat_int row, j, k, start, end;
    double *outp, *ap;
#ifdef USE_MPI
    int wrank;
    MPI_Comm mpicom = a->mpicom;
#ifdef USE_TASK
    const double *gather;
    MPI_Request sent = MPI_REQUEST_NULL, received = MPI_REQUEST_NULL;
#else // USE_TASK
    mat_int offset = 0;
    int rank, wsize;
    double *gather;
    MPI_Comm_size(mpicom, &wsize);
#endif // USE_TASK
#else // USE_MPI
    const double *gather;
    assert(lin == a->columns);
    assert(lout == a->rows);
#endif // USE_MPI
    TIME_TYPE t0, t1;
#ifdef USE_TASK
    TaskList *task;
    int taskCount = a->taskCount;
#else
    int taskCount = 1;
#endif
#ifdef USE_BLOCK
    const integer b1 = a->blockSize, b2=b1*b1;
    const doublereal one=1.0;
    const integer inc=1;
    const char trans='T';
#endif

    memset(out, 0, lout * sizeof(double));
    SET_TIME(t0);
#ifdef USE_MPI
    MPI_Comm_rank(mpicom, &wrank);
#endif
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
#else // USE_MPI
        gather = in;
#endif // USE_MPI
        ADD_TIME(a->mpiTime, t1, t0);

        /* For a given task, perform matrix-vector multiply */

#ifdef USE_TASK
        start = task->start;
        end = task->end;
#else // USE_TASK
        start = 0;
        end = a->iCount;
#endif // USE_TASK

        /* Loop over rows

          Including the loop over tasks inside the omp parallel
          block crashed. */
#pragma omp parallel default(shared) private(row, outp, j, k, ap)
        {
            /* Using schedule(static, 1) and disabling dynamic teams
               significantly degraded performance.  */
#pragma omp for
            for(row=start; row<end; row++) {
                outp = out + a->i[row];
#ifdef USE_BLOCK
                for(k=a->ip[row], ap=a->value + b2*a->ip[row];
                    k<a->ip[row+1]; k++, ap += b2) {
                    j = a->j[k];
                    DGEMV(&trans, &b1, &b1, &one,
                          ap, &b1, gather + j, &inc, &one,
                          outp, &inc);
                }
#else
                for(k=a->ip[row], ap=a->value; k<a->ip[row+1]; k++) {
                    j = a->j[k];
                    *outp += ap[k] * gather[j];
                }
#endif
            } // loop over rows
        } // omp parallel
        ADD_TIME(a->localTime, t0, t1);
#ifdef USE_MPI
#ifdef USE_TASK
        /* Verify that any current messages are complete
           before starting the next task or exiting.
	   Subsequent code may modify the 'in' array. */
        if(task->sendTo != -1)
            MPI_Wait(&sent, MPI_STATUS_IGNORE);
        if(task->receiveFrom != -1)
            MPI_Wait(&received, MPI_STATUS_IGNORE);
#endif // USE_TASK
#endif // USE_MPI
    } // loop over tasks

    ADD_TIME(a->mpiTime, t1, t0);
#endif  // MKL or not
    a->count += 1;
}
