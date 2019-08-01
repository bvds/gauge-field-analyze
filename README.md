# Chroma-generated gauge configuration analysis

## Running

```
gauge-analyze [input_config] [output]
 - [input_config]  the QDP config stored in the SciDAC file format (.lime)
 - [output] name of the output file in Mathematica format
```

Converts a configuration stored in the SciDAC file format (also known as
[lime](https://github.com/usqcd-software/c-lime)) to Mathematica-readable
format.

The input file is currently assumed to be created by the `purgaug` program in
[chroma](https://github.com/JeffersonLab/chroma). This is because of the XML
path to the string containing the lattice sizes is of a certain format. This
will be extended in the future. The program only runs in serial, and
running it in parallel will cause the program to abort.

## Installation

Install [QDP++](https://github.com/usqcd-software/qdpxx)
and packages `autoconf` and `g++`.

In the root directory, generate the configuration files:

    autoreconf -f

To get configuration options:

    ./configure --help

To build:

    ./congfigure [options]
    make


## Bibliography

 - [`writeopenqcd.cc`](https://rqcd.ur.de:8443/regensburg-lattice/chroma/blob/master/lib/io/writeopenqcd.cc)
 - [qdp-to-openqcd](https://github.com/Irubataru/qdp-to-openqcd)
 - Sara Collins, Regensburg:
   - [Usage example for `purgage`](https://homepages.uni-regensburg.de/~cos14742/lqcd-1/exercise5/extras/purgaug.html)
   - [slides](https://homepages.uni-regensburg.de/~cos14742/lqcd-1/exercise5/extras/slides.pdf)

# MATLAB linear system solvers #

To find stationary points of the action for some lattice configuration
we need a Krylov space method for solving symmtric-indefinite
singular (or incompatible) linear systems.

The MINRES-QLP algorithm can handle this case, however the associated [MATLAB
code](https://www.mathworks.com/matlabcentral/fileexchange/42419-minres-qlp) contains some errors. Also, I cleaned up the code a bit:

* Initialization error in `minresxxxM` (lines 746, 747)

* Condition `flag != flag0` (line 543) means that MINRES is never called.

* Replace `length(x)>0` with `~isempty(x)`

* Value of w_{k-3} (lines 557 to 559) is never used.

* `gamal3` is never used.

* Initial value is overwritten:  `gamal2`, `u`, and `wl2`.

* Added explicit `end` to each function.

Modified versions:  [`matlab-lib/minresqlpFull.m`](matlab-lib/minresqlpFull.m)
and [`matlab-lib/minresqlp0.m`](matlab-lib/minresqlp0.m) (no MINRES version).
I used this mainly to validate my export of the code to Mathematica.
