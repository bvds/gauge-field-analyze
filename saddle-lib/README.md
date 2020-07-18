#  Find the saddle point of a quadratic function

Use Krylov-space methods to find the saddle-point of a
quadratic function given the Hessian matrix, gradient vector,
and some number of linear constraints.  In addition, remove
eigenvectors where there is no nearby saddle-point.

In the case of lattice gauge theory, gauge transforms represent
directions where the action is invariant.  Since we want to find
points that are nearby in link-field space, we disallow directions
associated with gauge transforms.

1. Use MINRES/MINRES-QLP to project out shifts associated with infinitesimal
gauge transforms.
2. Use TRLan to identify shifts that are too large, orthogonalizing
with respect to infinitesimal gauge transforms.  That is, identify
shifts where the quadratic expansion is no longer valid.
3. Then use MINRES/MINRES-QLP to to find the shift, orthogonalizing with
respect to infinitesimal gauge transforms and the large shifts.

Since, in both cases, we need to modify the code to add
orthgonalization against previous vectors, we modify the original
routines, rather than use a version embedded in a library.
Eventually, one might combine steps 2 and 3 so that the Lanzos
tri-diagonalization is not repeated.

## Install

Install cJSON from GitHub in `/usr/local/lib`.
Then run `sudo ldconfig` to update the linker cache.

Download source code:
* <https://github.com/bvds/minres-qlp> is a customized version
  of [MINRES-QLP](https://web.stanford.edu/group/SOL/software/minresqlp/)
  from Stanford U.
* <https://github.com/bvds/nutrlan> is a customized version
  [nuTRLan-0.2](https://codeforge.lbl.gov/projects/trlan/) from LBL.

After downloading these, you will need to compile the
associated libraries.
Then modify `Makefile` to point to these libraries.

### BLAS libraries

The code has been tested against various versions of the
BLAS libraries.  It appears that OpenBLAS is the fastest.

1. OpenBLAS.  Compile from source (GitHub) and install in `/opt` (default).

2. Blis.  Compile from source (GitHub) and install in `/usr/local` (default).

3. Intel MKL. Install the libaray in `/opt` (default).
One can also link to the Intel Math Kernel Library
for the sparse matrix-vector multiplications.

The dynamic linker needs to find the BLAS library at runtime.
For example:
     *  Create file `/etc/ld.so.conf.d/mkl.conf`
        containing `/opt/intel/mkl/lib/intel64`
     *  run `sudo ldconfig`

The results are consistent with [an AMD Threadripper benchmark](
https://github.com/xianyi/OpenBLAS/issues/1461#issuecomment-469252560).


### LibRSB

I tried [`librsb`](http://librsb.sourceforge.net), but it
did not perform very well.  In addition, it cannot take advantage
of the color block structure of the matrices.

## Run

See the main program file `shifts.c` for instructions on
running the excutables `shifts` and `pshifts`.
