#  Apply Krylov-space methods to saddle-point step

1. Use MINRES/MINRES-QLP to project out shifts associated with infinitesimal
gauge transforms.
2. Use TRLan to identify shifts that are too large, orthogonalizing
with respect to infinitesimal gauge transforms.  That is, identify
shifts where the quadratic expansion is no longer valid.
3. Then use MINRES/MINRES-QLP to to find the shift, orthogonalizing with
respect to infinitesimal gauge transforms and the large shifts.

Since, in both cases, we need to modify the code to add
orthgonalization against previous vectors, we use the original
routines, rather than something embedded in a library.
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

### CentOs

On CentOs, had to provide a link for the LAPACK and BLAS libraries like this:  `sudo ln -s liblapack.so.3 liblapack.so`

### Intel MKL

One can link to the Intel Math Kernel Library
for the sparse matrix-vector multiplications.
*  Onstall the library in the default directory (`/opt`)
*  The dynamic linker needs to find the library:
   *  Create file `/etc/ld.so.conf.d/mkl.conf`
      containing `/opt/intel/mkl/lib/intel64`
   *  run `sudo ldconfig`

### LibRSB

I tried [`librsb`](http://librsb.sourceforge.net), but it
did not perform very well.  In addition, it cannot take advantage
of the color block structure of the matrices.
