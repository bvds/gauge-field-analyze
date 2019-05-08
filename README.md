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
   