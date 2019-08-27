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
I had trouble with the runtime finding library in `/usr/local/lib` and
ran `sudo ldconfig` to update the linker cache.

Download source code:

* [MINRES](https://web.stanford.edu/group/SOL/software/minres/) and
  [MINRES-QLP](https://web.stanford.edu/group/SOL/software/minresqlp/)
  from Stanford U.

* [nutrlan](https://codeforge.lbl.gov/projects/trlan/) from LBL.

and Make the libraries.  Modify `Makefile` to point to these libraries.
