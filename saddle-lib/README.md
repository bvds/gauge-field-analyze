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
* <https://github.com/bvds/minres-qlp> is a customized version
  of [MINRES-QLP](https://web.stanford.edu/group/SOL/software/minresqlp/)
  from Stanford U.
* <https://github.com/bvds/nutrlan> is a customized version
  [nuTRLan-0.2](https://codeforge.lbl.gov/projects/trlan/) from LBL.

After downloading these, you will need to compile the
associated libraries.

* In `minres-qlp`, specify the user-suppied
orthogonalizer ...

* For `nutrlan`, you will need to specify
the name of the user-supplied orthogonalizer
by adding the flag `-DUSER_ORTHO=dynamicProject` to `OPT`
in `Make.inc`.

Finally, modify `Makefile` to point to these libraries.

