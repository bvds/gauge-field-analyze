# SU(N) gauge field analysis

This repository contains code to analyze SU(N) lattice gauge
theory field configurations.  The strategy is to
look at 3 and 4 dimensions in the large N limit, paying special
attention to saddle points of the action.  It has three components:

1.  A number of Mathematica notebooks (`gauge.nb` being the
central file) containing test code and a number of analyses.
2.  A Mathematica library [`mma-lib`](mma-lib) for analyzing
gauge field configurations.
3.  C++ code for converting Chroma-generated gauge
configurations into Mathematica format.  See "gauge-analyze" below.
4.  C code for finding the saddle-point of a quadratic function,
given the Hessian matrix, the gradient vector, and some number of linear
constraints.  See [`saddle-lib`](saddle-lib).

## Chroma

The official release of [chroma](https://github.com/JeffersonLab/chroma)
does not properly handle large N or 3 spacetime dimensions.
Instead use my [fork of chroma](https://github.com/bvds/chroma).
Also, this fork includes an extension `TRANS_GAUGEBC` that sets
fixed boundary conditions in two directions.

## gauge-analyze

### Install

Install [QDP++](https://github.com/usqcd-software/qdpxx)
and packages `autoconf` and `g++`.

In the root directory, generate the configuration files:

    autoreconf -f
    automake --add-missing 

To get configuration options:

    ./configure --help

In particular, note the `--with-qdp` option.  To build, one
can create a separate build directory like this:

    mkdir 3-3-build/ ; cd 3-3-build
    ../configure --with-qdp=/usr/local/qdp++/3-3/
    make


### Run

```
ganalyze [input_config] [output]
 - [input_config]  the QDP config stored in the SciDAC file format (.lime)
 - [output] name of the output file in Mathematica format
```

Converts a configuration stored in the SciDAC file format (also known as
[lime](https://github.com/usqcd-software/c-lime)) to Mathematica-readable
format.

The input file is generated by the `purgaug` program from `chroma`
as discussed above.  The program only runs for a single thread, non-parallel.

Example commands to generate gauge configurations and create an
associated Mathematica input file:

    cd ../data/3-3/
    # create file periodic-16-28-in.xml
    ../../chroma/3-3-build/mainprogs/main/purgaug -i periodic-16-28-in.xml -o periodic-16-28-out.xml
    ../../gauge-field-analyze/3-3-build/programs/ganalyze periodic-16-28.ini.xml5 periodic-16-28-5.m

Or, if there is a lot of them:

    perl -e 'for (1..100) {system "../../gauge-field-analyze/programs/ganalyze periodic-16-28.ini.xml$_ periodic-16-28-$_.m"}'


## Bibliography

 - [`writeopenqcd.cc`](https://rqcd.ur.de:8443/regensburg-lattice/chroma/blob/master/lib/io/writeopenqcd.cc)
 - [qdp-to-openqcd](https://github.com/Irubataru/qdp-to-openqcd)
 - Sara Collins, Regensburg:
   - [Usage example for `purgaug`](https://homepages.uni-regensburg.de/~cos14742/lqcd-1/exercise5/extras/purgaug.html)
   - [slides](https://homepages.uni-regensburg.de/~cos14742/lqcd-1/exercise5/extras/slides.pdf)

# MATLAB linear system solvers #

To find stationary points of the action for some lattice configuration
we need a Krylov space method for solving symmtric-indefinite
singular (or incompatible) linear systems.

The MINRES-QLP algorithm can handle this case, however the associated [MATLAB
code](https://www.mathworks.com/matlabcentral/fileexchange/42419-minres-qlp) contains some errors. Also, I cleaned up the code a bit:

* Condition `flag != flag0` (line 543) means that MINRES is never called.

* Initialization error in `minresxxxM` (lines 746, 747)

* Replace `length(x)>0` with `~isempty(x)`

* Value of w_{k-3} (lines 557 to 559) is never used.

* `gamal3` is never used.

* Initial value is overwritten:  `gamal2`, `u`, and `wl2`.

* Added explicit `end` to each function.

MATLAB versions:

* [`matlab-lib/minresqlp.m`](matlab-lib/minresqlp.m) Starts with MINRES
  with switch to MINRES-QLP.

* [`matlab-lib/minresqlp0.m`](matlab-lib/minresqlp0.m) Uses MINRES-QLP only.

* [`matlab-lib/minres1.m`](matlab-lib/minres1.m) Uses MINRES only.  Should
 be functionally equivalent to the [original MINRES function](https://web.stanford.edu/group/SOL/software/minres/).

Mathematica versions:

* [`mma-lib/minresqlp.m`](mma-lib/minresqlp.m) Starts with MINRES
  with switch to MINRES-QLP.

* [`mma-lib/minresqlp0.m`](mma-lib/minresqlp0.m) Uses MINRES-QLP only.

* [`mma-lib/minres.m`](mma-lib/minres.m) Port of the [original MINRES function](https://web.stanford.edu/group/SOL/software/minres/).

* [`mma-lib/minres1.m`](mma-lib/minres1.m) Version of MINRES-QLP with MINRES only.  Should be equivalent to [original MINRES function](https://web.stanford.edu/group/SOL/software/minres/).
