#
#   Stand-alone program, maybe called from Mathematica
#
MINRESDIR = ../../linear-algebra/minres-qlp
TRLANDIR = ../../linear-algebra/nutrlan
MKLROOT = /opt/intel/mkl
PRIMMEDIR = ../../linear-algebra/primme
#
#
#                                 -g -Wall -O3
CFLAGS = -fPIC -I$(TRLANDIR)/INC -g -Wall -Wextra -pedantic -std=gnu11 -fopenmp
LFLAGS =  -L$(TRLANDIR) -L/usr/local/lib -lcjson -lgomp
SLFLAGS = -lnutrlan
#
#  Need to add fortran libraries needed by MINRES-QLP
#  Use "ldd minresqlptest_mpi" to view libraries.
PLFLAGS = -lmpi_mpifh -lmpi_usempif08 -lnutrlan_mpi
#
#  On AMD Epyc, add flags suggested by
#  http://www.prace-ri.eu/best-practice-guide-amd-epyc
#
#  Possible additional flags, some associated with prefetching:
#   -ffast-math -ftree-vectorize -fprefetch-loop-arrays -mprefetchwt1
#  These need to be tested emprically because they might
#  interfere with hardware prefetching.
#
#  glibc 2.22+ a vectorized version of the libm library is available.
#  Use -lmvec instead of -lm, also need flags -ftree-vectorize
#  -funsafe-math-optimizations and -ffast-math.
#
ifeq ($(shell lscpu|grep -m 1 -c -i epyc),1)
CFLAGS += -O3 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer
LLAPACK = -lflame
else
LLAPACK = -llapack
endif
#
#   Use block-sparse coordinate format for sparse matrices
#
CFLAGS += -DUSE_BLOCK
#
#   Use task manager to control sparse matrix-vector products
#
CFLAGS += -DUSE_TASK
#
#        BLAS choices
#
ifeq (1,1) # Use OpenBlas
CFLAGS += -I/opt/OpenBLAS/include/ -DUSE_OPENBLAS
LFLAGS += $(LLAPACK) -L/opt/OpenBLAS/lib/ -lopenblas
#
else ifeq (1,0) #   Use Intel MKL
#   Also, use MKL for sparse matrix storage if the flags
#   USE_MPI, USE_BLOCK, and USE_TASK have not been set.
#   This gives a substantial gain relative to default BLAS and
#   plain matrix storage.
#
#  Flags generated by tool:
#  cd /opt/intel/mkl/tools; ./mkl_link_tool -interactive
#
#  compiler flags:
CFLAGS += -DUSE_MKL -m64 -I$(MKLROOT)/include
#
#  Linker flags:
#  Looks like dynamic linking doesn't work well if the routine
#  is called by Mathematica.
#  Use static links:
LFLAGS +=  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -ldl
#  Dynamic:
#LFLAGS += -L$(MKLROOT)/lib/intel64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -Wl,--no-as-needed -lmkl_rt -lgomp -ldl
#LFLAGS += -L$(MKLROOT)/lib/intel64  -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -ldl
#
else ifeq (1,0)  # Use Blis
CFLAGS += -DUSE_BLIS
LFLAGS += $(LLAPACK) -lblis
#
else
LFLAGS += $(LLAPACK) -lblas
endif
# These need to be at the end
LFLAGS += -lgfortran -lm
#
#
#  Optionally, use Jacobi-Davidson library PRIMME
#
#  I abandoned this attempt because adding a function user_ortho
#  to enforce the gauge transform constraint appears to be difficult.
ifeq (1, 0)
CFLAGS += -DUSE_PRIMME -I$(PRIMMEDIR)/include
LFLAGS += -L$(PRIMMEDIR)/lib -lprimme
endif
#

MINRESOBJS = $(MINRESDIR)/minresqlpDataModule.o  $(MINRESDIR)/minresqlpBlasModule.o
SMINRESOBJS = $(MINRESOBJS) $(MINRESDIR)/minresqlpModule.o 
PMINRESOBJS = $(MINRESOBJS) $(MINRESDIR)/minresqlpModule_mpi.o $(MINRESDIR)/parallelBlasModule_mpi.o

SOBJS = shifts.o dynamic.o eigensystem.o cutoff.o linear.o mmio.o sort.o matrix.o
POBJS = shifts_mpi.o dynamic_mpi.o eigensystem_mpi.o cutoff_mpi.o
POBJS += linear_mpi.o mmio.o sort.o matrix_mpi.o

MPICC = mpicc
MPIFC = mpifc

shifts.h: fortran.h

shifts: $(SOBJS)
	$(CC) $(CFLAGS) -o $@ $(SOBJS) $(SMINRESOBJS) $(SLFLAGS) $(LFLAGS)

pshifts: $(POBJS)
	$(MPICC) $(CFLAGS) -o $@ $(POBJS) $(PMINRESOBJS) $(PLFLAGS) $(LFLAGS)

%.o: %.c shifts.h mmio.h Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

%_mpi.o: %.c shifts.h Makefile
	$(MPICC) -c -DUSE_MPI -o $@ $< $(CFLAGS)

prun:
	mpirun -np 2 xterm -e gdb -ex run --args ./pshifts hess-grad-gauge.json junk.out

clean:
	rm shifts *.o
