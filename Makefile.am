#
# Unlike chroma, using a single top-level Makefile
#
#  BvdS:  temporary to get stuff running
#
OPENQCDDIR=/home/bvds/lattice/openQCD-1.6

#
# All the include files (for distribution)
#
include_HEADERS = \
     lib/openqcd_archive.hpp    lib/qdp_converters.hpp \
     lib/openqcd_utilities.hpp  lib/qdp_utilities.hpp

#
# Include path and any other compiler flags 
#
AM_CXXFLAGS  = -I$(srcdir)/lib -I$(top_builddir)/lib -I$(OPENQCDDIR)/include
AM_CXXFLAGS += @CXXFLAGS@ @QDPXX_CXXFLAGS@

AM_CFLAGS = -I$(srcdir)/lib -I$(top_builddir)/lib 
AM_CFLAGS += @CFLAGS@

lib_LIBRARIES = libganalyze.a

libganalyze_a_SOURCES = \
    lib/openqcd_archive.cpp    lib/qdp_converters.cpp \
    lib/openqcd_utilities.cpp  lib/qdp_utilities.cpp 

#
# Include Path
#
INCFLAGS=-I$(srcdir)/programs -I$(srcdir)/lib -I$(top_builddir)/lib
#
AM_LDFLAGS = -L$(BUILDDIR)/lib  @QDPXX_LDFLAGS@ @LDFLAGS@
LDADD = -lganalyze @QDPXX_LIBS@ @LIBS@ 

bin_PROGRAMS = programs/ganalyze
programs_ganalyze_SOURCES = programs/gauge-analyze.cpp
